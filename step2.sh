cd /mnt/sdb/transcriptome/elakatophyceae_pep/pp/deredundancy/0212/cross/OrthoFinder/

# 1) 先从 OG 列表提取 fasta
while read -r og; do
  cp "Orthogroup_Sequences/${og}.fa" SOG468_fa/
done < Orthogroups/SOG_5of5_1copy.txt

# 2) 确认确实提取到了文件
ls SOG468_fa/*.fa | wc -l

# 3) 并行 mafft + hmmbuild
JOBS=196
find SOG468_fa -name "*.fa" | xargs -I{} -P ${JOBS} bash -c '
f="$1"
b=$(basename "$f" .fa)
mafft --thread 1 --auto "$f" > "SOG468_aln/${b}.aln.fa" 2> "logs/${b}.mafft.log" &&
hmmbuild "SOG468_hmm/${b}.hmm" "SOG468_aln/${b}.aln.fa" > "logs/${b}.hmmbuild.log" 2>&1
' _ {}

# 4) 合并并压库
ls SOG468_hmm/*.hmm | wc -l
cat SOG468_hmm/*.hmm > SOG468.hmmdb
hmmpress SOG468.hmmdb
