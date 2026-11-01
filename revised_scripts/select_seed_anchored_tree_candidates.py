#!/usr/bin/env python3
"""Select one putative ortholog per OTU from multi-hit HMM gene trees.

HMM scores are used only to retain a bounded candidate pool upstream.  For each
non-seed OTU, this script compares *all* retained candidates in the ML gene
tree and selects the copy with the smallest median patristic distance to the
five fixed seed orthologs.  A score rank is used exclusively to break an exact
tree-distance tie.  The outputs provide a record for every selection so that
the procedure can be independently audited.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path
from statistics import mean, median

from Bio import Phylo


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-audit", required=True, type=Path)
    parser.add_argument("--tree-dir", required=True, type=Path)
    parser.add_argument("--fasta-dir", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--n-taxa", required=True, type=int)
    parser.add_argument("--min-taxa", required=True, type=int)
    parser.add_argument("--seed", action="append", required=True)
    return parser.parse_args()


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    header = None
    sequence: list[str] = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records[header] = "".join(sequence)
                header, sequence = line[1:].split()[0], []
            else:
                sequence.append(line)
    if header is not None:
        records[header] = "".join(sequence)
    return records


def otu(header: str) -> str:
    if "|" not in header:
        raise ValueError(f"Header lacks an OTU delimiter: {header}")
    return header.split("|", 1)[0]


def numeric_rank(row: dict[str, str]) -> int:
    value = row.get("candidate_rank_after_quality_filter", "")
    try:
        return int(value)
    except ValueError:
        return 0


def numeric_score(row: dict[str, str]) -> float:
    try:
        return float(row.get("full_sequence_bitscore", "nan"))
    except ValueError:
        return float("nan")


def main() -> None:
    args = parse_args()
    seeds = set(args.seed)
    if len(seeds) < 2:
        raise SystemExit("At least two fixed seed OTUs are required")

    metadata: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)
    with args.candidate_audit.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            metadata[row["locus"]][row["tree_header"]] = row

    out_dir = args.out_dir
    if out_dir.exists() and any(out_dir.iterdir()):
        raise SystemExit(f"Refusing to overwrite non-empty output directory: {out_dir}")
    raw_dir = out_dir / "raw_selected_ge_threshold"
    raw_dir.mkdir(parents=True, exist_ok=True)

    sequence_rows: list[dict[str, object]] = []
    locus_rows: list[dict[str, object]] = []
    selected_loci: list[str] = []

    for locus in sorted(metadata):
        tree_path = args.tree_dir / f"{locus}.treefile"
        fasta_path = args.fasta_dir / f"{locus}.fas"
        if not tree_path.is_file() or not fasta_path.is_file():
            locus_rows.append({
                "locus": locus, "status": "missing_tree_or_fasta", "n_tree_leaves": 0,
                "n_otus_input": 0, "n_selected_otus": 0, "coverage": 0,
                "n_fixed_seed_otus": 0, "n_multicopy_nonseed_otus": 0,
                "n_tree_guided_changes_from_rank1": 0, "mean_assignment_margin": "",
                "median_assignment_margin": "", "eligible_for_realign": 0,
            })
            continue

        tree = Phylo.read(str(tree_path), "newick")
        terminals = {leaf.name: leaf for leaf in tree.get_terminals()}
        candidates_by_otu: dict[str, list[str]] = defaultdict(list)
        for header in terminals:
            if header not in metadata[locus]:
                continue
            candidates_by_otu[otu(header)].append(header)

        present_seeds = seeds & set(candidates_by_otu)
        seed_headers = []
        for seed in sorted(present_seeds):
            fixed = [h for h in candidates_by_otu[seed] if "|SEED_" in h]
            if len(fixed) == 1:
                seed_headers.append(fixed[0])

        selected: dict[str, tuple[str, float, float | None, int]] = {}
        if len(seed_headers) == len(seeds):
            for seed_header in seed_headers:
                selected[otu(seed_header)] = (seed_header, 0.0, None, 1)
            for sample, headers in sorted(candidates_by_otu.items()):
                if sample in seeds:
                    continue
                scored = []
                for header in headers:
                    distances = [tree.distance(terminals[header], terminals[s]) for s in seed_headers]
                    scored.append((median(distances), mean(distances), header))
                # Topological/branch-length proximity is primary.  HMM rank and
                # bitscore are deterministic tie-breakers only.
                scored.sort(key=lambda x: (
                    x[0], x[1], numeric_rank(metadata[locus][x[2]]) or 10**9,
                    -numeric_score(metadata[locus][x[2]]), x[2],
                ))
                best = scored[0]
                runner_up = scored[1] if len(scored) > 1 else None
                margin = None
                if runner_up is not None:
                    margin = (runner_up[0] - best[0]) / max(best[0], 1e-12)
                selected[sample] = (best[2], best[0], margin, len(scored))

        fasta = read_fasta(fasta_path)
        n_selected = len(selected)
        coverage = n_selected / args.n_taxa
        eligible = int(
            len(seed_headers) == len(seeds)
            and n_selected >= args.min_taxa
            and all(header in fasta for header, _, _, _ in selected.values())
        )
        margins = [m for _, _, m, n in selected.values() if m is not None and n > 1]
        multicopy_otus = sum(len(headers) > 1 for sample, headers in candidates_by_otu.items() if sample not in seeds)
        changes = 0
        for sample, (header, _, _, _) in selected.items():
            if sample in seeds:
                continue
            if numeric_rank(metadata[locus][header]) != 1:
                changes += 1

        locus_rows.append({
            "locus": locus,
            "status": "evaluated" if len(seed_headers) == len(seeds) else "missing_fixed_seed",
            "n_tree_leaves": len(terminals),
            "n_otus_input": len(candidates_by_otu),
            "n_selected_otus": n_selected,
            "coverage": f"{coverage:.6f}",
            "n_fixed_seed_otus": len(seed_headers),
            "n_multicopy_nonseed_otus": multicopy_otus,
            "n_tree_guided_changes_from_rank1": changes,
            "mean_assignment_margin": f"{mean(margins):.6f}" if margins else "",
            "median_assignment_margin": f"{median(margins):.6f}" if margins else "",
            "eligible_for_realign": eligible,
        })

        for sample, (header, distance, margin, candidate_count) in sorted(selected.items()):
            row = metadata[locus][header]
            sequence_rows.append({
                "locus": locus,
                "otu": sample,
                "selected_tree_header": header,
                "source": row["source"],
                "original_target": row["original_target"],
                "candidate_rank_after_quality_filter": row["candidate_rank_after_quality_filter"],
                "full_sequence_bitscore": row["full_sequence_bitscore"],
                "union_hmm_coverage": row["union_hmm_coverage"],
                "n_candidates_compared_in_otu": candidate_count,
                "median_patristic_distance_to_five_fixed_seeds": f"{distance:.6f}",
                "relative_margin_to_second_tree_candidate": "" if margin is None else f"{margin:.6f}",
                "selected_by_tree_not_hmm_rank1": int(sample not in seeds and numeric_rank(row) != 1),
                "eligible_locus": eligible,
            })

        if eligible:
            with (raw_dir / f"{locus}.fas").open("w") as handle:
                for sample, (header, _, _, _) in sorted(selected.items()):
                    handle.write(f">{sample}\n{fasta[header]}\n")
            selected_loci.append(locus)

    if not locus_rows:
        raise SystemExit("No loci were evaluated")
    out_dir.mkdir(parents=True, exist_ok=True)
    with (out_dir / "locus_selection_audit.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(locus_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(locus_rows)
    with (out_dir / "sequence_selection_audit.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(sequence_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(sequence_rows)
    with (out_dir / "selected_loci_pre_realign.txt").open("w") as handle:
        handle.write("\n".join(selected_loci) + ("\n" if selected_loci else ""))
    with (out_dir / "selection_parameters.tsv").open("w") as handle:
        handle.write("parameter\tvalue\n")
        handle.write(f"n_taxa\t{args.n_taxa}\n")
        handle.write(f"minimum_taxa_pre_realign\t{args.min_taxa}\n")
        handle.write(f"fixed_seed_otus\t{','.join(sorted(seeds))}\n")
        handle.write("selection_rule\tFor each non-seed OTU, all retained HMM candidates were compared in the ML gene tree. The candidate with the smallest median patristic distance to the five fixed seed orthologs was selected; HMM rank and bitscore were used only to resolve exact tree-distance ties.\n")

    print(f"Loci evaluated: {len(locus_rows)}")
    print(f"Loci meeting pre-realignment coverage threshold: {len(selected_loci)}")
    print(f"Selected raw FASTA directory: {raw_dir}")


if __name__ == "__main__":
    main()
