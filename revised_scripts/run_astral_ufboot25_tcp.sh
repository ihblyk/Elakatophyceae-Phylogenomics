#!/usr/bin/env bash
set -Eeuo pipefail

RUN="/mnt/sdb/Revised.ele/05_68taxa_seedanchored_rerun_excluding_SRR9052848_v1"
SEL="${RUN}/05_seed_anchored_selection"
TREE_DIR="${SEL}/05_final_gene_trees_MFP"
OUT="${SEL}/06_astral_ufboot25_tcp_constraint"
STATUS_DIR="${RUN}/status"

EXPECTED_LOCI=388
PYTHON="/home/ps/local/share/mamba/envs/ortho/bin/python"
PREP="/mnt/sdb/Revised.ele/04_multihit_orthology_screen/scripts/prepare_astral_collapsed_gene_trees.py"
ASTRAL="/home/ps/ASTRAL/Astral/astral.5.7.8.jar"
CONSTRAINED_ASTRAL="/mnt/sdb/Revised.ele/tools/astral_constrained_5.6.9/Astral/astral.5.6.9.jar"

die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 1
}

for executable in "$PYTHON" java; do
  command -v "$executable" >/dev/null 2>&1 || die "Missing executable: $executable"
done

for required in "$PREP" "$ASTRAL" "$CONSTRAINED_ASTRAL"; do
  [[ -s "$required" ]] || die "Missing required file: $required"
done

mkdir -p "$STATUS_DIR"

# Do not consume a partial gene-tree collection while the parallel IQ-TREE stage is active.
while :; do
  n_trees=$(find "$TREE_DIR" -maxdepth 1 -type f -name '*.treefile' | wc -l)
  if [[ "$n_trees" -eq "$EXPECTED_LOCI" ]] && [[ -f "${STATUS_DIR}/05_final_gene_trees.complete" ]]; then
    break
  fi
  printf '[%s] waiting for final gene trees: %s/%s\n' "$(date '+%F %T')" "$n_trees" "$EXPECTED_LOCI"
  sleep 60
done

[[ ! -e "$OUT" ]] || die "Refusing to overwrite existing output directory: $OUT"
mkdir -p "$OUT/collapsed_gene_trees"

COMBINED="${OUT}/gene_trees_ufboot25.tre"
AUDIT="${OUT}/ufboot25_collapse_audit.tsv"
CONSTRAINT="${OUT}/TCP_monophyly_constraint.tre"

printf '%s\n' '(SRR12875398_cdhit,SRR16849199_cdhit,SRR14413272_cdhit);' > "$CONSTRAINT"

"$PYTHON" "$PREP" \
  --tree-dir "$TREE_DIR" \
  --out-dir "${OUT}/collapsed_gene_trees" \
  --out-combined "$COMBINED" \
  --out-audit "$AUDIT" \
  --threshold 25

# The inherited helper labels outputs as collapsed35; rename them to the actual threshold.
find "${OUT}/collapsed_gene_trees" -maxdepth 1 -type f -name '*.collapsed35.tree' -print0 | while IFS= read -r -d '' tree; do
  mv "$tree" "${tree%.collapsed35.tree}.collapsed25.tree"
done

n_collapsed=$(find "${OUT}/collapsed_gene_trees" -maxdepth 1 -type f -name '*.collapsed25.tree' | wc -l)
[[ "$n_collapsed" -eq "$EXPECTED_LOCI" ]] || die "Expected ${EXPECTED_LOCI} collapsed gene trees, found ${n_collapsed}"
[[ "$(grep -c ';' "$COMBINED")" -eq "$EXPECTED_LOCI" ]] || die "Combined gene-tree count differs from ${EXPECTED_LOCI}"

java -Xmx96G -jar "$ASTRAL" \
  -i "$COMBINED" \
  -o "${OUT}/astral_ufboot25_unconstrained_v5.7.8.tree" \
  -t 3 \
  > "${OUT}/astral_ufboot25_unconstrained_v5.7.8.stdout.log" \
  2> "${OUT}/astral_ufboot25_unconstrained_v5.7.8.stderr.log"

# This matched v5.6.9 run makes the score comparison with the constrained search version-controlled.
java -Xmx96G -jar "$CONSTRAINED_ASTRAL" \
  -i "$COMBINED" \
  -o "${OUT}/astral_ufboot25_unconstrained_v5.6.9.tree" \
  -t 3 \
  > "${OUT}/astral_ufboot25_unconstrained_v5.6.9.stdout.log" \
  2> "${OUT}/astral_ufboot25_unconstrained_v5.6.9.stderr.log"

java -Xmx96G -jar "$CONSTRAINED_ASTRAL" \
  -i "$COMBINED" \
  -j "$CONSTRAINT" \
  -o "${OUT}/astral_ufboot25_TCP_constrained_v5.6.9.tree" \
  -t 3 \
  > "${OUT}/astral_ufboot25_TCP_constrained_v5.6.9.stdout.log" \
  2> "${OUT}/astral_ufboot25_TCP_constrained_v5.6.9.stderr.log"

# Score both species-tree hypotheses with the same current ASTRAL implementation.
for tree in \
  "${OUT}/astral_ufboot25_unconstrained_v5.7.8.tree" \
  "${OUT}/astral_ufboot25_TCP_constrained_v5.6.9.tree"; do
  label=$(basename "$tree" .tree)
  java -Xmx96G -jar "$ASTRAL" \
    -i "$COMBINED" \
    -q "$tree" \
    -t 3 \
    -o "${OUT}/${label}.quartet_scored.tree" \
    > "${OUT}/${label}.quartet_score.stdout.log" \
    2> "${OUT}/${label}.quartet_score.stderr.log"
done

"$PYTHON" - "$OUT" <<'PY'
from pathlib import Path
from Bio import Phylo
import sys

out = Path(sys.argv[1])
tcp = {
    "SRR12875398_cdhit",  # Pithophora roettleri
    "SRR16849199_cdhit",  # Trentepohlia annulata
    "SRR14413272_cdhit",  # Chlorocladiella medogensis
}
trees = sorted(out.glob("astral_ufboot25_*.tree"))

with (out / "TCP_constraint_validation.tsv").open("w") as handle:
    handle.write("tree\ttcp_exact_edge_bipartition\n")
    for path in trees:
        tree = Phylo.read(path, "newick")
        all_taxa = {tip.name for tip in tree.get_terminals()}
        observed = False
        for clade in tree.find_clades(order="preorder"):
            side = {tip.name for tip in clade.get_terminals()}
            if side == tcp or all_taxa - side == tcp:
                observed = True
                break
        handle.write(f"{path.name}\t{str(observed).lower()}\n")
        if path.name.endswith("TCP_constrained_v5.6.9.tree") and not observed:
            raise SystemExit("TCP constraint was not recovered as an exact edge bipartition")
PY

{
  printf 'parameter\tvalue\n'
  printf 'taxa\t68\n'
  printf 'final_loci\t%s\n' "$EXPECTED_LOCI"
  printf 'ufboot_collapse_rule\t<25\n'
  printf 'constraint\tPithophora_roettleri+Trentepohlia_annulata+Chlorocladiella_medogensis monophyly\n'
  printf 'unconstrained_astral\t5.7.8\n'
  printf 'constraint_aware_astral\t5.6.9 (-j)\n'
  printf 'constraint_newick\t%s\n' "$(cat "$CONSTRAINT")"
} > "${OUT}/analysis_parameters.tsv"

touch "${STATUS_DIR}/06_astral_ufboot25_tcp.complete"
printf 'Completed ASTRAL analyses: %s\n' "$OUT"
