#!/usr/bin/env bash
# Re-align and quality-filter the tree-guided, one-sequence-per-OTU selection.
set -Eeuo pipefail

RUN="${1:?Usage: run_realign_treeguided_selection.sh <selection-output-directory>}"
INPUT_DIR="$RUN/raw_selected_ge_threshold"
ALN_DIR="$RUN/02_mafft"
TRIM_DIR="$RUN/03_trimal"
FINAL_DIR="$RUN/04_final_ge80_len100"
LOG_DIR="$RUN/logs_realign"
STATS="$RUN/realign_filter_audit.tsv"

JOBS="${JOBS:-100}"
MAFFT_THREADS="${MAFFT_THREADS:-2}"
N_TAXA="${N_TAXA:-69}"
MIN_TAXA="${MIN_TAXA:-56}"
MIN_ALN_LEN="${MIN_ALN_LEN:-100}"

MAFFT="${MAFFT:-/home/ps/local/share/mamba/envs/ortho/bin/mafft}"
TRIMAL="${TRIMAL:-/home/ps/local/share/mamba/envs/ortho/bin/trimal}"

[[ -d "$INPUT_DIR" ]] || { echo "ERROR: missing input directory: $INPUT_DIR" >&2; exit 1; }
[[ -x "$MAFFT" ]] || { echo "ERROR: MAFFT not executable: $MAFFT" >&2; exit 1; }
[[ -x "$TRIMAL" ]] || { echo "ERROR: trimAl not executable: $TRIMAL" >&2; exit 1; }

for d in "$ALN_DIR" "$TRIM_DIR" "$FINAL_DIR" "$LOG_DIR"; do
  if [[ -e "$d" ]] && find "$d" -mindepth 1 -print -quit | grep -q .; then
    echo "ERROR: refusing to overwrite non-empty directory: $d" >&2
    exit 1
  fi
  mkdir -p "$d"
done

mapfile -d '' FILES < <(find "$INPUT_DIR" -maxdepth 1 -type f -name '*.fas' -print0 | sort -z)
[[ ${#FILES[@]} -gt 0 ]] || { echo "ERROR: no FASTA files found in $INPUT_DIR" >&2; exit 1; }

export ALN_DIR TRIM_DIR LOG_DIR MAFFT TRIMAL MAFFT_THREADS
run_one() {
  local input="$1"
  local locus
  locus=$(basename "$input" .fas)
  "$MAFFT" --thread "$MAFFT_THREADS" --auto "$input" \
    > "$ALN_DIR/${locus}.fas" 2> "$LOG_DIR/${locus}.mafft.log"
  "$TRIMAL" -in "$ALN_DIR/${locus}.fas" -out "$TRIM_DIR/${locus}.fas" \
    -gt 0.75 -resoverlap 0.75 -seqoverlap 70 \
    > "$LOG_DIR/${locus}.trimal.log" 2>&1
}
export -f run_one

printf 'Input loci: %s; parallel jobs: %s; MAFFT threads/job: %s\n' "${#FILES[@]}" "$JOBS" "$MAFFT_THREADS"
printf '%s\0' "${FILES[@]}" | xargs -0 -r -n 1 -P "$JOBS" bash -c 'run_one "$1"' _

python3 - "$TRIM_DIR" "$FINAL_DIR" "$STATS" "$N_TAXA" "$MIN_TAXA" "$MIN_ALN_LEN" <<'PY'
import csv
import glob
import os
import shutil
import sys

trim_dir, final_dir, stats_file, n_total, min_taxa, min_len = sys.argv[1:]
n_total, min_taxa, min_len = map(int, (n_total, min_taxa, min_len))

def read_fasta(path):
    records = []
    header, chunks = None, []
    with open(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    records.append((header, ''.join(chunks)))
                header, chunks = line[1:].split()[0], []
            else:
                chunks.append(line)
    if header is not None:
        records.append((header, ''.join(chunks)))
    return records

rows = []
kept = []
for path in sorted(glob.glob(os.path.join(trim_dir, '*.fas'))):
    records = read_fasta(path)
    headers = [h for h, _ in records]
    duplicate_headers = len(headers) != len(set(headers))
    lengths = {len(s) for _, s in records}
    aln_len = next(iter(lengths)) if len(lengths) == 1 else -1
    nseq = len(records)
    min_ungapped = min((len(s.replace('-', '')) for _, s in records), default=0)
    median_ungapped = sorted((len(s.replace('-', '')) for _, s in records))
    median_ungapped = median_ungapped[len(median_ungapped)//2] if median_ungapped else 0
    coverage = nseq / n_total
    keep = int((not duplicate_headers) and nseq >= min_taxa and aln_len >= min_len)
    rows.append({
        'locus': os.path.basename(path).removesuffix('.fas'),
        'n_otus_after_trimal': nseq,
        'coverage_after_trimal': f'{coverage:.6f}',
        'alignment_length': aln_len,
        'min_ungapped_sequence_length': min_ungapped,
        'median_ungapped_sequence_length': median_ungapped,
        'duplicate_headers': int(duplicate_headers),
        'kept_ge80_len100': keep,
    })
    if keep:
        shutil.copy2(path, os.path.join(final_dir, os.path.basename(path)))
        kept.append(os.path.basename(path).removesuffix('.fas'))

with open(stats_file, 'w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
    writer.writeheader()
    writer.writerows(rows)
with open(os.path.join(os.path.dirname(stats_file), 'final_loci_ge80_len100.txt'), 'w') as handle:
    handle.write('\n'.join(kept) + ('\n' if kept else ''))

print('Trimmed loci:', len(rows))
print('Final loci (>=%.0f%% occupancy; alignment >=%d aa):' % (100*min_taxa/n_total, min_len), len(kept))
print('Final directory:', final_dir)
PY

echo "Trimmed loci: $(find "$TRIM_DIR" -maxdepth 1 -type f -name '*.fas' | wc -l)"
echo "Final loci: $(find "$FINAL_DIR" -maxdepth 1 -type f -name '*.fas' | wc -l)"
