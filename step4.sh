cd /mnt/sdb/transcriptome/elakatophyceae_pep/pp/deredundancy/0212/cross/Orthofinder
mkdir -p step4_tables step4_loci_fa

python - <<'PY'
import os, glob, math
from collections import defaultdict

dom_dir = "hmmsearch_domtbl"
pep_dir = "/mnt/sdb/transcriptome/elakatophyceae_pep/pp/deredundancy/0212/cross/all"   # 你的pep文件目录
out_dir = "step4_loci_fa"
os.makedirs(out_dir, exist_ok=True)

# 1) 解析domtbl：每个 sample × query 只保留最佳命中（最小iE，若并列取更高score）
best = {}  # (sample, query) -> (target, iE, score)
for fp in glob.glob(os.path.join(dom_dir, "*.domtbl")):
    sample = os.path.basename(fp).replace(".domtbl","")
    with open(fp, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            x = line.split()
            target = x[0]
            query  = x[3]
            iE     = float(x[12])  # i-Evalue
            score  = float(x[13])  # domain score
            k = (sample, query)
            if (k not in best) or (iE < best[k][1]) or (math.isclose(iE, best[k][1]) and score > best[k][2]):
                best[k] = (target, iE, score)

# 保存best-hit表
with open("step4_tables/best_hits.tsv", "w") as o:
    o.write("sample\tquery\ttarget\tiEvalue\tscore\n")
    for (s,q),(t,e,sc) in sorted(best.items()):
        o.write(f"{s}\t{q}\t{t}\t{e}\t{sc}\n")

# 2) 统计每个sample需要提取哪些蛋白ID
need = defaultdict(set)
for (s,q),(t,e,sc) in best.items():
    need[s].add(t)

# 3) 从pep中提取所需序列
seqs = {}  # (sample, id) -> seq
for pep in glob.glob(os.path.join(pep_dir, "*.pep")):
    sample = os.path.basename(pep).replace(".pep","")
    if sample not in need:
        continue
    wanted = need[sample]
    hid, buf = None, []
    with open(pep, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                if hid is not None and hid in wanted:
                    seqs[(sample, hid)] = "".join(buf)
                hid = line[1:].strip().split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if hid is not None and hid in wanted:
            seqs[(sample, hid)] = "".join(buf)

# 4) 按query输出per-locus fasta
missing = 0
counts = defaultdict(int)
for (s,q),(t,e,sc) in sorted(best.items()):
    qsafe = q.replace("/", "_").replace(" ", "_")
    fa = os.path.join(out_dir, f"{qsafe}.fa")
    seq = seqs.get((s,t))
    if not seq:
        missing += 1
        continue
    with open(fa, "a") as o:
        o.write(f">{s}|{t}\n{seq}\n")
    counts[qsafe] += 1

with open("step4_tables/per_locus_counts.tsv", "w") as o:
    o.write("locus\tnseq\n")
    for k,v in sorted(counts.items()):
        o.write(f"{k}\t{v}\n")

print("best-hit pairs:", len(best))
print("loci generated:", len(counts))
print("missing seq pairs:", missing)
print("best table: step4_tables/best_hits.tsv")
print("counts: step4_tables/per_locus_counts.tsv")
PY
