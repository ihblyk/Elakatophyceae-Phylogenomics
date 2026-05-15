#!/bin/bash
set -euo pipefail

# =========================
# Config
# =========================
WORKDIR="cd /mnt/sdb/transcriptome/elakatophyceae_pep/pp/deredundancy/0212/cross/OrthoFinder"
IN_DIR="${WORKDIR}/step6_trim_renamed"

TREE_DIR="${WORKDIR}/step7_iqtree"
TS_IN="${WORKDIR}/step7_treeshrink_in"
TS_OUT="${WORKDIR}/step7_treeshrink_out"
TS_FA="${WORKDIR}/step7_treeshrink_fa"

FINAL_DIR="${WORKDIR}/step7_loci_final_ge80_len100"
STATS_TSV="${WORKDIR}/step7_loci_filter_stats.tsv"
LOG_DIR="${WORKDIR}/logs_step7"

JOBS=220            # 并行任务数（按机器调）
IQ_THREADS=1       # 每个IQ-TREE任务线程
IQ_MODEL_MODE="MFP_RESTRICTED" # 可改 MFP，但更慢
TS_Q=0.05          # TreeShrink参数
MIN_OCC=0.80       # 覆盖度阈值
MIN_LEN=100        # 长度阈值（aa）


# =========================
# Check
# =========================
command -v iqtree2 >/dev/null 2>&1 || { echo "ERROR: iqtree2 not found"; exit 1; }
command -v run_treeshrink.py >/dev/null 2>&1 || { echo "ERROR: run_treeshrink.py not found"; exit 1; }
command -v python3 >/dev/null 2>&1 || { echo "ERROR: python3 not found"; exit 1; }

mkdir -p "$TREE_DIR" "$TS_IN" "$TS_OUT" "$TS_FA" "$FINAL_DIR" "$LOG_DIR"
rm -f "${TREE_DIR}"/* "${TS_FA}"/* "${FINAL_DIR}"/* 2>/dev/null || true
rm -rf "${TS_IN:?}/"* "${TS_OUT:?}/"* 2>/dev/null || true
export TREE_DIR LOG_DIR IQ_THREADS JOBS

# =========================
# 1) IQ-TREE gene trees
# =========================
find "$IN_DIR" -maxdepth 1 -type f -name "*.fa" -print0 | \
xargs -0 -n1 -P ${JOBS} bash -c '
f="$1"
b=$(basename "$f" .fa)
iqtree2 -s "$f" \
  -m MFP \
  -mset LG,WAG,Q.plant \
  -bb 1000 \
  -T ${IQ_THREADS} -quiet -redo \
  -pre "${TREE_DIR}/${b}" > "${LOG_DIR}/${b}.iqtree.log" 2>&1
' _

# =========================
# 2) Prepare TreeShrink input
# =========================
for fa in "$IN_DIR"/*.fa; do
  b=$(basename "$fa" .fa)
  tree="${TREE_DIR}/${b}.treefile"
  [[ -s "$tree" ]] || { echo "WARN: missing treefile for $b"; continue; }
  mkdir -p "${TS_IN}/${b}"
  cp "$fa"   "${TS_IN}/${b}/input.fasta"
  cp "$tree" "${TS_IN}/${b}/input.tree"
done

# =========================
# 3) TreeShrink
# =========================
run_treeshrink.py \
  -i "$TS_IN" \
  -t input.tree \
  -a input.fasta \
  -q "$TS_Q" \
  -o "$TS_OUT" \
  > "${LOG_DIR}/treeshrink.log" 2>&1

# =========================
# 4) Collect TreeShrink FASTA
# =========================
# 兼容不同版本输出位置
count=0
while IFS= read -r -d '' f; do
  b=$(basename "$(dirname "$f")")
  cp "$f" "${TS_FA}/${b}.fa"
  count=$((count+1))
done < <(find "$TS_OUT" -type f -name "output.fasta" -print0)

if [[ $count -eq 0 ]]; then
  while IFS= read -r -d '' f; do
    b=$(basename "$(dirname "$f")")
    cp "$f" "${TS_FA}/${b}.fa"
    count=$((count+1))
  done < <(find "$TS_IN" -type f -name "output.fasta" -print0)
fi

echo "TreeShrink FASTA collected: $count"

# =========================
# 5) Filter loci by occupancy and length
#    occupancy >= 0.80, length >= 100 aa
# =========================
python3 - <<PY
import os, glob

in_dir = "${TS_FA}"
final_dir = "${FINAL_DIR}"
stats_tsv = "${STATS_TSV}"
min_occ = float("${MIN_OCC}")
min_len = int("${MIN_LEN}")

os.makedirs(final_dir, exist_ok=True)

def read_fasta(fp):
    recs = []
    h = None
    buf = []
    with open(fp) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if h is not None:
                    recs.append((h, "".join(buf)))
                h = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if h is not None:
            recs.append((h, "".join(buf)))
    return recs

# 用原始输入目录估计总taxa数（期望覆盖分母）
all_taxa = set()
for fp in glob.glob(os.path.join("${IN_DIR}", "*.fa")):
    for h, s in read_fasta(fp):
        all_taxa.add(h)
n_total = len(all_taxa)
if n_total == 0:
    raise SystemExit("No taxa found in input alignments.")

kept = 0
rows = []

for fp in sorted(glob.glob(os.path.join(in_dir, "*.fa"))):
    locus = os.path.basename(fp)
    recs = read_fasta(fp)

    # 去除重复taxon（保留第一条）
    seen = set()
    uniq = []
    for h, s in recs:
        if h in seen:
            continue
        seen.add(h)
        uniq.append((h, s))

    nseq = len(uniq)
    aln_len = len(uniq[0][1]) if nseq else 0
    occ = nseq / n_total
    ok = (occ >= min_occ and aln_len >= min_len and nseq > 0)

    rows.append((locus, nseq, aln_len, occ, int(ok)))

    if ok:
        out = os.path.join(final_dir, locus)
        with open(out, "w") as w:
            for h, s in uniq:
                w.write(f">{h}\n{s}\n")
        kept += 1

with open(stats_tsv, "w") as o:
    o.write("locus\tnseq\taln_len\toccupancy\tkept\n")
    for r in rows:
        o.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]:.4f}\t{r[4]}\n")

print("Total taxa:", n_total)
print("TreeShrink loci:", len(rows))
print("Kept loci:", kept)
print("Final dir:", final_dir)
print("Stats:", stats_tsv)
PY

echo "Input loci:  $(find "$IN_DIR" -maxdepth 1 -name '*.fa' | wc -l)"
echo "Final loci:  $(find "$FINAL_DIR" -maxdepth 1 -name '*.fa' | wc -l)"
echo "Done."
