cd /mnt/sdb/transcriptome/elakatophyceae_pep/pp/deredundancy/0212/cross/OrthoFinder
mkdir -p step6_aln step6_trim logs_step6

JOBS=196
find step5_loci_ge80 -name "*.fa" | xargs -I{} -P ${JOBS} bash -c '
f="$1"; b=$(basename "$f" .fa)
mafft --thread 1 --auto "$f" > "step6_aln/${b}.aln.fa" 2> "logs_step6/${b}.mafft.log" &&
trimal -in "step6_aln/${b}.aln.fa" -out "step6_trim/${b}.trim.fa" -gt 0.70 -cons 50 > "logs_step6/${b}.trimal.log" 2>&1
' _ {}

ls step6_aln/*.aln.fa | wc -l
ls step6_trim/*.trim.fa | wc -l

