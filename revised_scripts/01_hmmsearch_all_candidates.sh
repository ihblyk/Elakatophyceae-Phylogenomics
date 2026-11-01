#!/usr/bin/env bash
set -Eeuo pipefail

ROOT=/mnt/sdb/Revised.ele
OUT="$ROOT/04_multihit_orthology_screen"
PEP_DIR="$ROOT/01_inputs/phylogenomics_69_pep"
DB="$ROOT/03_active_hmm_library_5of5_1copy/SOG_5of5_1copy_478.hmmdb"
HMMSEARCH=/home/ps/local/share/mamba/envs/rna_denovo/bin/hmmsearch
TBL="$OUT/hmmsearch_tbl"
DOM="$OUT/hmmsearch_domtbl"
LOG="$OUT/logs/hmmsearch"
DONE="$OUT/run/hmmsearch_done"
JOBS=100
THREADS=2
EVAL=1e-5
DOMEVAL=1e-5

mkdir -p "$TBL" "$DOM" "$LOG" "$DONE"

run_one() {
  local q="$1"
  local base
  base=$(basename "$q" .pep)
  "$HMMSEARCH" --cpu "$THREADS" --noali --notextw \
    -E "$EVAL" --domE "$DOMEVAL" \
    --tblout "$TBL/$base.tbl" \
    --domtblout "$DOM/$base.domtbl" \
    "$DB" "$q" > "$LOG/$base.hmmsearch.log" 2>&1
  touch "$DONE/$base.done"
}
export HMMSEARCH DB TBL DOM LOG DONE THREADS EVAL DOMEVAL
export -f run_one

find "$PEP_DIR" -maxdepth 1 -type f -name '*.pep' -print0 | sort -z | \
  xargs -0 -r -n 1 -P "$JOBS" bash -c 'run_one "$1"' _

input_count=$(find "$PEP_DIR" -maxdepth 1 -type f -name '*.pep' | wc -l)
done_count=$(find "$DONE" -maxdepth 1 -type f -name '*.done' | wc -l)
[ "$input_count" -eq "$done_count" ] || {
  echo "ERROR: completed $done_count of $input_count searches" >&2
  exit 1
}
printf 'completed_samples\t%s\n' "$done_count" > "$OUT/manifest/hmmsearch_completion.tsv"
printf 'candidate_threshold\tfull_sequence_E<=%s; domain_E<=%s\n' "$EVAL" "$DOMEVAL" >> "$OUT/manifest/hmmsearch_completion.tsv"
printf 'parallel_slots\t%s\nthreads_per_task\t%s\n' "$JOBS" "$THREADS" >> "$OUT/manifest/hmmsearch_completion.tsv"
