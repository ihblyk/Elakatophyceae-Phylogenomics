#!/usr/bin/env python3
"""Collapse weak UFBoot branches in IQ-TREE consensus gene trees for ASTRAL."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from Bio import Phylo


def ufboot(node_name: str | None) -> float | None:
    """Read the right-hand (UFBoot) value from IQ-TREE SH-aLRT/UFBoot labels."""
    if not node_name:
        return None
    candidate = node_name.split("/")[-1]
    try:
        return float(candidate)
    except ValueError:
        return None


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree-dir", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--out-combined", required=True, type=Path)
    parser.add_argument("--out-audit", required=True, type=Path)
    parser.add_argument("--threshold", type=float, default=35.0)
    args = parser.parse_args()

    # IQ-TREE's .treefile retains SH-aLRT/UFBoot labels, whereas .contree
    # contains only the UFBoot consensus values.
    source_files = sorted(args.tree_dir.glob("*.treefile"))
    if not source_files:
        raise SystemExit(f"No IQ-TREE .treefile files in {args.tree_dir}")
    if args.out_dir.exists() and any(args.out_dir.iterdir()):
        raise SystemExit(f"Refusing to overwrite non-empty directory: {args.out_dir}")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.out_combined.parent.mkdir(parents=True, exist_ok=True)
    args.out_audit.parent.mkdir(parents=True, exist_ok=True)

    audit_rows = []
    combined_trees = []
    for source in source_files:
        tree = Phylo.read(str(source), "newick")
        weak = []
        labeled = 0
        for clade in tree.find_clades(order="postorder"):
            if clade is tree.root or clade.is_terminal():
                continue
            value = ufboot(clade.name)
            if value is not None:
                labeled += 1
                if value < args.threshold:
                    weak.append(clade)
        # Postorder permits nested low-support nodes to be collapsed safely.
        collapsed = 0
        for clade in weak:
            try:
                tree.collapse(clade)
                collapsed += 1
            except ValueError:
                # A parent may already have been collapsed; this has the same
                # intended outcome for the nested child.
                pass
        # ASTRAL needs clean Newick topology.  The support labels have already
        # been used for collapse and are retained in the original IQ-TREE files
        # plus the audit table, so remove residual SH-aLRT/UFBoot labels here.
        for clade in tree.find_clades():
            if not clade.is_terminal():
                clade.name = None
        destination = args.out_dir / source.name.replace(".treefile", ".collapsed35.tree")
        Phylo.write(tree, str(destination), "newick")
        combined_trees.append(tree)
        audit_rows.append({
            "locus": source.stem,
            "n_tips": len(tree.get_terminals()),
            "n_labeled_internal_branches": labeled,
            "n_branches_collapsed_below_threshold": collapsed,
            "threshold_ufboot": args.threshold,
        })

    Phylo.write(combined_trees, str(args.out_combined), "newick")
    with args.out_audit.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(audit_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(audit_rows)
    print(f"Processed gene trees: {len(combined_trees)}")
    print(f"Combined input: {args.out_combined}")
    print(f"Audit: {args.out_audit}")


if __name__ == "__main__":
    main()
