#!/usr/bin/env python3
"""Build and select seed-anchored, tree-pruned orthology candidates.

This script deliberately retains multiple HMM candidates per sample and locus.
It therefore separates permissive candidate retrieval from the subsequent
tree-based, single-copy orthology decision made by PhyloPyPruner.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path
from typing import Iterable


SEED_TO_OTU = {
    "Charabraunii": "Charabraunii_cdhit",
    "Chlorokybusatmophyticus": "Chlorokybusatmophyticus_cdhit",
    "ERR3487426": "ERR3487426_cdhit",
    "Mesostigmaviride": "Mesostigmaviride_cdhit",
    "SRR5949008": "SRR5949008_cdhit",
}
SEED_OTUS = frozenset(SEED_TO_OTU.values())


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    header: str | None = None
    sequence: list[str] = []

    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records[header] = "".join(sequence).replace("*", "")
                header = line[1:].split()[0]
                sequence = []
            else:
                sequence.append(line)

    if header is not None:
        records[header] = "".join(sequence).replace("*", "")
    return records


def write_fasta(path: Path, records: Iterable[tuple[str, str]]) -> None:
    with path.open("w") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n{sequence}\n")


def safe_identifier(value: str) -> str:
    # Keep exactly one delimiter in output headers for PhyloPyPruner OTU parsing.
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")


def as_float(value: str) -> float:
    return 0.0 if value == "0" else float(value)


def read_loci(path: Path) -> list[str]:
    loci = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if len(loci) != len(set(loci)):
        raise SystemExit("The seed locus list contains duplicate IDs.")
    if len(loci) != 478:
        raise SystemExit(f"Expected 478 strict seed loci, found {len(loci)}.")
    return loci


def get_seed_records(seed_path: Path) -> dict[str, tuple[str, str, str]]:
    raw_records = read_fasta(seed_path)
    if len(raw_records) != 5:
        raise SystemExit(f"{seed_path} has {len(raw_records)} records; expected five.")

    records: dict[str, tuple[str, str, str]] = {}
    for raw_header, sequence in raw_records.items():
        seed_name, sep, identifier = raw_header.partition("@")
        if not sep or seed_name not in SEED_TO_OTU:
            raise SystemExit(f"Unrecognized seed header: {raw_header}")
        otu = SEED_TO_OTU[seed_name]
        records[otu] = (raw_header, sequence, identifier)

    if set(records) != SEED_OTUS:
        raise SystemExit(f"Seed mapping mismatch in {seed_path.name}.")
    return records


def build(args: argparse.Namespace) -> None:
    screen_root = Path(args.screen_root)
    pep_dir = Path(args.pep_dir)
    seed_dir = Path(args.seed_dir)
    candidate_table = Path(args.candidate_table)
    locus_list = Path(args.locus_list)
    output_dir = Path(args.output_dir)
    input_dir = output_dir / "fasta"

    if output_dir.exists():
        raise SystemExit(f"Output directory already exists: {output_dir}")
    input_dir.mkdir(parents=True)

    loci = read_loci(locus_list)
    locus_set = set(loci)

    raw_count: dict[str, int] = defaultdict(int)
    quality: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    with candidate_table.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            locus = row["locus"]
            if locus not in locus_set:
                continue
            raw_count[locus] += 1
            coverage = float(row["union_hmm_coverage"])
            evalue = as_float(row["full_sequence_evalue"])
            bitscore = float(row["full_sequence_bitscore"])
            if coverage < args.min_hmm_coverage or evalue > args.max_evalue:
                continue
            row["_coverage"] = coverage
            row["_evalue"] = evalue
            row["_bitscore"] = bitscore
            quality[(row["sample"], locus)].append(row)

    candidate_audit: list[list[object]] = []
    locus_stats: list[list[object]] = []
    retained_total = 0
    # Each PEP file is reused by many loci. Cache it once to avoid repeatedly
    # parsing the same transcriptome-scale file during candidate assembly.
    pep_cache: dict[str, dict[str, str]] = {}
    for locus in loci:
        seeds = get_seed_records(seed_dir / f"{locus}.fa")
        by_sample: dict[str, list[dict[str, object]]] = {}
        quality_count = 0
        relative_count = 0
        for (sample, candidate_locus), rows in quality.items():
            if candidate_locus != locus:
                continue
            # The five seed OTUs already have a verified strict-SOG sequence.
            # Do not let additional HMM matches from those reference proteomes
            # compete with that known ortholog during tree-based pruning.
            if sample in SEED_OTUS:
                continue
            quality_count += len(rows)
            best = max(float(row["_bitscore"]) for row in rows)
            retained = [
                row
                for row in rows
                if float(row["_bitscore"]) >= args.min_relative_bitscore * best
            ]
            retained.sort(
                key=lambda row: (
                    -float(row["_bitscore"]),
                    float(row["_evalue"]),
                    str(row["target"]),
                )
            )
            relative_count += len(retained)
            by_sample[sample] = retained[: args.max_candidates]

        records: list[tuple[str, str]] = []
        seen_headers: set[str] = set()
        selected_candidate_count = 0
        candidate_otus: set[str] = set()
        for sample in sorted(by_sample):
            pep_path = pep_dir / f"{sample}.pep"
            if not pep_path.is_file():
                raise SystemExit(f"Missing active PEP file: {pep_path}")
            if sample not in pep_cache:
                pep_cache[sample] = read_fasta(pep_path)
            pep_records = pep_cache[sample]
            for rank, row in enumerate(by_sample[sample], start=1):
                target = str(row["target"])
                if target not in pep_records:
                    raise SystemExit(f"{target} from {pep_path.name} is absent from its PEP file.")
                sequence = pep_records[target]
                # The fixed seed is retained once, rather than being duplicated as an
                # identical HMM candidate in the same OTU.
                if sample in seeds and sequence == seeds[sample][1]:
                    continue
                header = f"{sample}|CAND_{safe_identifier(target)}"
                if header in seen_headers:
                    continue
                seen_headers.add(header)
                records.append((header, sequence))
                selected_candidate_count += 1
                candidate_otus.add(sample)
                candidate_audit.append(
                    [
                        locus,
                        header,
                        "hmm_candidate",
                        sample,
                        target,
                        rank,
                        len(sequence),
                        row["hmm_length"],
                        row["full_sequence_evalue"],
                        row["full_sequence_bitscore"],
                        row["union_hmm_coverage"],
                        row["reported_domains"],
                    ]
                )

        for otu in sorted(SEED_OTUS):
            raw_header, sequence, identifier = seeds[otu]
            header = f"{otu}|SEED_{safe_identifier(identifier)}"
            if header in seen_headers:
                raise SystemExit(f"Duplicate fixed-seed header in {locus}: {header}")
            seen_headers.add(header)
            records.append((header, sequence))
            candidate_audit.append(
                [
                    locus,
                    header,
                    "fixed_seed_reference",
                    otu,
                    raw_header,
                    0,
                    len(sequence),
                    "",
                    "",
                    "",
                    "",
                    "",
                ]
            )

        write_fasta(input_dir / f"{locus}.fas", records)
        retained_total += selected_candidate_count
        locus_stats.append(
            [
                locus,
                raw_count[locus],
                quality_count,
                relative_count,
                selected_candidate_count,
                len(candidate_otus),
                len(records),
                5,
            ]
        )

    with (output_dir / "candidate_sequence_audit.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "locus",
                "tree_header",
                "source",
                "sample",
                "original_target",
                "candidate_rank_after_quality_filter",
                "sequence_length",
                "hmm_length",
                "full_sequence_evalue",
                "full_sequence_bitscore",
                "union_hmm_coverage",
                "reported_domains",
            ]
        )
        writer.writerows(candidate_audit)

    with (output_dir / "per_locus_candidate_input_stats.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "locus",
                "raw_domain_candidates",
                "quality_candidates",
                "relative_score_candidates_before_top_n",
                "retained_hmm_candidates_after_top_n_and_seed_deduplication",
                "active_taxa_with_retained_hmm_candidate",
                "total_tree_records_including_fixed_seeds",
                "fixed_seed_references",
            ]
        )
        writer.writerows(locus_stats)

    with (output_dir / "parameters.tsv").open("w") as handle:
        handle.write("parameter\tvalue\n")
        handle.write(f"seed_loci\t{len(loci)}\n")
        handle.write(f"minimum_union_hmm_coverage\t{args.min_hmm_coverage}\n")
        handle.write(f"maximum_full_sequence_evalue\t{args.max_evalue}\n")
        handle.write(
            "minimum_relative_bitscore_among_quality_candidates\t"
            f"{args.min_relative_bitscore}\n"
        )
        handle.write(f"maximum_candidates_per_sample_locus\t{args.max_candidates}\n")
        handle.write(
            "seed_otu_handling\tEach seed OTU is represented only by its original "
            "strict single-copy reference sequence; additional HMM candidates from "
            "seed proteomes are excluded before tree inference.\n"
        )
        handle.write(f"retained_hmm_candidates\t{retained_total}\n")
        handle.write(
            "selection_note\tAll qualifying candidates were retained up to the top-N cap; "
            "no single best HMM hit was selected at this stage.\n"
        )

    print(f"Loci written: {len(loci)}")
    print(f"Retained HMM candidates: {retained_total}")
    print(f"FASTA directory: {input_dir}")


def otu_from_header(header: str) -> str:
    if "|" in header:
        return header.split("|", 1)[0]
    return header.split("@", 1)[0]


def select(args: argparse.Namespace) -> None:
    pp_output = Path(args.pp_output)
    output_dir = Path(args.output_dir)
    locus_list = Path(args.locus_list)
    loci = read_loci(locus_list)
    groups = defaultdict(list)
    for path in sorted(pp_output.glob("*.fas")):
        locus, separator, _ = path.name.partition("_pruned_")
        if not separator:
            continue
        groups[locus].append(path)

    if output_dir.exists():
        if any(output_dir.iterdir()):
            raise SystemExit(f"Output directory is not empty: {output_dir}")
    else:
        output_dir.mkdir(parents=True)

    summary_rows: list[list[object]] = []
    selected_sequence_rows: list[list[object]] = []
    selected_loci: list[str] = []
    for locus in loci:
        evaluated: list[dict[str, object]] = []
        for path in groups[locus]:
            records = read_fasta(path)
            headers = list(records)
            otus = [otu_from_header(header) for header in headers]
            alignment_lengths = {len(sequence) for sequence in records.values()}
            alignment_length = next(iter(alignment_lengths)) if len(alignment_lengths) == 1 else 0
            evaluated.append(
                {
                    "path": path,
                    "records": records,
                    "otus": otus,
                    "unique": len(otus) == len(set(otus)),
                    "seed_count": len(set(otus) & SEED_OTUS),
                    "alignment_length": alignment_length,
                    "min_ungapped_length": min(
                        (len(sequence.replace("-", "")) for sequence in records.values()),
                        default=0,
                    ),
                }
            )

        qualified = [
            item
            for item in evaluated
            if item["unique"]
            and SEED_OTUS.issubset(set(item["otus"]))
            and len(item["otus"]) >= args.min_taxa
            and item["alignment_length"] >= args.min_alignment_length
        ]

        if qualified:
            qualified.sort(
                key=lambda item: (
                    -len(item["otus"]),
                    -int(item["alignment_length"]),
                    str(item["path"]),
                )
            )
            selected = qualified[0]
            output_path = output_dir / f"{locus}.fas"
            normalized = [(otu_from_header(header), sequence) for header, sequence in selected["records"].items()]
            write_fasta(output_path, normalized)
            selected_loci.append(locus)
            for header, sequence in selected["records"].items():
                selected_sequence_rows.append(
                    [locus, otu_from_header(header), header, len(sequence.replace("-", ""))]
                )
            summary_rows.append(
                [
                    locus,
                    "selected",
                    selected["path"].name,
                    len(selected["otus"]),
                    selected["alignment_length"],
                    selected["min_ungapped_length"],
                    selected["seed_count"],
                    len(evaluated),
                    "",
                ]
            )
            continue

        anchor_groups = [item for item in evaluated if SEED_OTUS.issubset(set(item["otus"]))]
        if not evaluated:
            reason = "no_pruned_group"
        elif not anchor_groups:
            reason = "no_group_retained_all_five_seed_otus"
        elif not any(item["unique"] for item in anchor_groups):
            reason = "seed_anchored_group_not_single_copy"
        elif not any(len(item["otus"]) >= args.min_taxa for item in anchor_groups):
            reason = "seed_anchored_group_below_taxon_coverage"
        else:
            reason = "seed_anchored_group_below_alignment_length"
        best_anchor = max(anchor_groups, key=lambda item: len(item["otus"]), default=None)
        summary_rows.append(
            [
                locus,
                "discarded",
                best_anchor["path"].name if best_anchor else "",
                len(best_anchor["otus"]) if best_anchor else 0,
                best_anchor["alignment_length"] if best_anchor else 0,
                best_anchor["min_ungapped_length"] if best_anchor else 0,
                best_anchor["seed_count"] if best_anchor else 0,
                len(evaluated),
                reason,
            ]
        )

    parent = output_dir.parent
    with (parent / "locus_selection_summary.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "locus",
                "status",
                "selected_or_best_anchor_group",
                "n_otus",
                "alignment_length",
                "minimum_ungapped_sequence_length",
                "seed_otus_present",
                "n_pruned_groups_evaluated",
                "reason_if_discarded",
            ]
        )
        writer.writerows(summary_rows)

    with (parent / "selected_sequence_manifest.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["locus", "final_otu", "source_tree_header", "ungapped_sequence_length"])
        writer.writerows(selected_sequence_rows)

    with (parent / "selected_loci.txt").open("w") as handle:
        for locus in selected_loci:
            handle.write(f"{locus}\n")

    print(f"Selected loci: {len(selected_loci)}")
    print(f"Final FASTA directory: {output_dir}")


def audit(args: argparse.Namespace) -> None:
    """Summarize every tree-pruned candidate group before final selection."""
    pp_output = Path(args.pp_output)
    output_path = Path(args.output)
    if output_path.exists():
        raise SystemExit(f"Audit output already exists: {output_path}")

    rows: list[list[object]] = []
    for path in sorted(pp_output.glob("*.fas")):
        locus, separator, _ = path.name.partition("_pruned_")
        if not separator:
            continue
        records = read_fasta(path)
        headers = list(records)
        otus = [otu_from_header(header) for header in headers]
        alignment_lengths = {len(sequence) for sequence in records.values()}
        seed_otus = sorted(set(otus) & SEED_OTUS)
        seed_reference_records = sum("|SEED_" in header for header in headers)
        rows.append(
            [
                locus,
                path.name,
                len(headers),
                len(set(otus)),
                int(len(headers) == len(set(otus))),
                len(seed_otus),
                ",".join(seed_otus),
                seed_reference_records,
                next(iter(alignment_lengths)) if len(alignment_lengths) == 1 else 0,
                min(
                    (len(sequence.replace("-", "")) for sequence in records.values()),
                    default=0,
                ),
            ]
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "locus",
                "pruned_group",
                "n_sequences",
                "n_otus",
                "is_single_copy_by_otu",
                "n_seed_otus",
                "seed_otus",
                "n_fixed_seed_records",
                "alignment_length",
                "minimum_ungapped_sequence_length",
            ]
        )
        writer.writerows(rows)
    print(f"Pruned groups audited: {len(rows)}")
    print(f"Audit table: {output_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    build_parser = subparsers.add_parser("build")
    build_parser.add_argument("--screen-root", required=True)
    build_parser.add_argument("--pep-dir", required=True)
    build_parser.add_argument("--seed-dir", required=True)
    build_parser.add_argument("--candidate-table", required=True)
    build_parser.add_argument("--locus-list", required=True)
    build_parser.add_argument("--output-dir", required=True)
    build_parser.add_argument("--min-hmm-coverage", type=float, default=0.60)
    build_parser.add_argument("--max-evalue", type=float, default=1e-10)
    build_parser.add_argument("--min-relative-bitscore", type=float, default=0.70)
    build_parser.add_argument("--max-candidates", type=int, default=5)
    build_parser.set_defaults(func=build)

    select_parser = subparsers.add_parser("select")
    select_parser.add_argument("--pp-output", required=True)
    select_parser.add_argument("--locus-list", required=True)
    select_parser.add_argument("--output-dir", required=True)
    select_parser.add_argument("--min-taxa", type=int, default=56)
    select_parser.add_argument("--min-alignment-length", type=int, default=100)
    select_parser.set_defaults(func=select)

    audit_parser = subparsers.add_parser("audit")
    audit_parser.add_argument("--pp-output", required=True)
    audit_parser.add_argument("--output", required=True)
    audit_parser.set_defaults(func=audit)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
