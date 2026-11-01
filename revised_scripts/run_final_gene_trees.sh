#!/usr/bin/env bash
# Infer support-annotated ML gene trees for the tree-guided single-copy matrix.
set -Eeuo pipefail

RUN="${1:?Usage: run_final_gene_trees.sh <selection-output-directory>}"
INPUT_DIR="$RUN/04_final_ge80_len100"
TREE_DIR="$RUN/05_final_gene_trees_MFP"
LOG_DIR="$RUN/logs_final_gene_trees"

JOBS="${JOBS:-100}"
IQ_THREADS="${IQ_THREADS:-2}"
IQTREE="${IQTREE:-/home/ps/local/share/mamba/envs/qtp18s_iqtree2_240/bin/iqtree2}"

[[ -d "$INPUT_DIR" ]] || { echo "ERROR: missing input directory: $INPUT_DIR" >&2; exit 1; }
[[ -x "$IQTREE" ]] || { echo "ERROR: IQ-TREE not executable: $IQTREE" >&2; exit 1; }
for d in "$TREE_DIR" "$LOG_DIR"; do
  if [[ -e "$d" ]] && find "$d" -mindepth 1 -print -quit | grep -q .; then
    echo "ERROR: refusing to overwrite non-empty directory: $d" >&2
    exit 1
  fi
  mkdir -p "$d"
done

mapfile -d '' FILES < <(find "$INPUT_DIR" -maxdepth 1 -type f -name '*.fas' -print0 | sort -z)
[[ ${#FILES[@]} -gt 0 ]] || { echo "ERROR: no final alignments in $INPUT_DIR" >&2; exit 1; }

export TREE_DIR LOG_DIR IQTREE IQ_THREADS
run_iqtree() {
  local input="$1"
  local locus
  locus=$(basename "$input" .fas)
  "$IQTREE" -s "$input" \
    -m MFP -mset LG,WAG,JTT,Q.pfam,Q.plant \
    -bb 1000 -alrt 1000 -bnni -T "$IQ_THREADS" -redo -quiet \
    -pre "$TREE_DIR/$locus" > "$LOG_DIR/${locus}.iqtree.log" 2>&1
}
export -f run_iqtree

printf 'Final loci: %s; parallel jobs: %s; IQ-TREE threads/job: %s\n' "${#FILES[@]}" "$JOBS" "$IQ_THREADS"
printf '%s\0' "${FILES[@]}" | xargs -0 -r -n 1 -P "$JOBS" bash -c 'run_iqtree "$1"' _

tree_count=$(find "$TREE_DIR" -maxdepth 1 -type f -name '*.treefile' | wc -l)
printf 'Completed gene trees: %s/%s\n' "$tree_count" "${#FILES[@]}"
[[ "$tree_count" -eq "${#FILES[@]}" ]] || {
  echo "ERROR: one or more gene-tree jobs failed; inspect $LOG_DIR" >&2
  exit 1
}
