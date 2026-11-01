#!/usr/bin/env python3
"""Build a validated supermatrix, partition file, and TCP-guided start tree."""

from __future__ import annotations

import argparse
from pathlib import Path

from Bio import Phylo, SeqIO


RHODOPHYTA_TIPS = {
    "SRR6023957_cdhit",
    "SRR10150610_cdhit",
    "SRR12352384_cdhit",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignments", required=True, type=Path)
    parser.add_argument("--constraint-tree", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--expected-loci", type=int, default=388)
    parser.add_argument("--expected-taxa", type=int, default=68)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    files = sorted(args.alignments.glob("*.fas"))
    if len(files) != args.expected_loci:
        raise SystemExit(f"Expected {args.expected_loci} alignments, found {len(files)}")

    tree = Phylo.read(args.constraint_tree, "newick")
    taxa = [tip.name for tip in tree.get_terminals()]
    if len(taxa) != args.expected_taxa or len(taxa) != len(set(taxa)):
        raise SystemExit(f"Constraint tree must contain {args.expected_taxa} unique tips")

    concatenated = {taxon: [] for taxon in taxa}
    partitions: list[tuple[str, int, int, int, int]] = []
    observed_taxa: set[str] = set()
    coordinate = 1

    for file_path in files:
        records = list(SeqIO.parse(file_path, "fasta"))
        records_by_id = {record.id: str(record.seq).upper() for record in records}
        if len(records_by_id) != len(records):
            raise SystemExit(f"Duplicate sequence IDs in {file_path.name}")
        unexpected = set(records_by_id) - set(taxa)
        if unexpected:
            raise SystemExit(f"Unexpected taxa in {file_path.name}: {sorted(unexpected)}")
        lengths = {len(sequence) for sequence in records_by_id.values()}
        if len(lengths) != 1:
            raise SystemExit(f"Non-aligned sequences in {file_path.name}")
        length = lengths.pop()
        if length == 0:
            raise SystemExit(f"Empty alignment in {file_path.name}")

        observed_taxa.update(records_by_id)
        for taxon in taxa:
            concatenated[taxon].append(records_by_id.get(taxon, "-" * length))

        start = coordinate
        end = coordinate + length - 1
        partitions.append((file_path.stem, start, end, length, len(records_by_id)))
        coordinate = end + 1

    if observed_taxa != set(taxa):
        raise SystemExit("Taxon set in alignments does not exactly match the constraint tree")

    matrix_length = coordinate - 1
    fasta_path = args.outdir / "supermatrix_68taxa_388loci.fas"
    phylip_path = args.outdir / "supermatrix_68taxa_388loci.phy"
    partition_path = args.outdir / "partitions_388loci.nex"
    locus_audit_path = args.outdir / "locus_partition_audit.tsv"
    taxa_path = args.outdir / "taxa_order.txt"
    start_tree_path = args.outdir / "TCP_start_tree_Rhodophyta_rooted_68taxa.nwk"

    with fasta_path.open("w", encoding="ascii") as handle:
        for taxon in taxa:
            handle.write(f">{taxon}\n{''.join(concatenated[taxon])}\n")

    with phylip_path.open("w", encoding="ascii") as handle:
        handle.write(f"{len(taxa)} {matrix_length}\n")
        for taxon in taxa:
            handle.write(f"{taxon} {''.join(concatenated[taxon])}\n")

    with partition_path.open("w", encoding="ascii") as handle:
        handle.write("#NEXUS\n\nbegin sets;\n")
        for locus, start, end, _, _ in partitions:
            handle.write(f"  charset {locus} = {start}-{end};\n")
        handle.write("end;\n")

    with locus_audit_path.open("w", encoding="ascii") as handle:
        handle.write("locus\tstart\tend\talignment_length\tn_observed_taxa\n")
        for row in partitions:
            handle.write("\t".join(map(str, row)) + "\n")

    taxa_path.write_text("\n".join(taxa) + "\n", encoding="ascii")

    outgroup_node = tree.common_ancestor(*RHODOPHYTA_TIPS)
    outgroup_tips = {tip.name for tip in outgroup_node.get_terminals()}
    if outgroup_tips != RHODOPHYTA_TIPS:
        raise SystemExit("Rhodophyta outgroup is not an exact monophyletic clade")
    edge = outgroup_node.branch_length or 1.0
    tree.root_with_outgroup(outgroup_node, outgroup_branch_length=edge / 2.0)
    root_children = tree.root.clades
    if not any({tip.name for tip in child.get_terminals()} == RHODOPHYTA_TIPS for child in root_children):
        raise SystemExit("Failed to place Rhodophyta as one root child")
    for clade in tree.find_clades():
        clade.branch_length = None
        if not clade.is_terminal():
            clade.confidence = None
            clade.name = None
    Phylo.write(tree, start_tree_path, "newick")

    print(f"loci={len(files)}")
    print(f"taxa={len(taxa)}")
    print(f"matrix_sites={matrix_length}")
    print(f"fasta={fasta_path}")
    print(f"phylip={phylip_path}")
    print(f"partitions={partition_path}")
    print(f"start_tree={start_tree_path}")


if __name__ == "__main__":
    main()
