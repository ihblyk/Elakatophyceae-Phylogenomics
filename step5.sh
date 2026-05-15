cd /mnt/sdb/transcriptome/elakatophyceae_pep/pp/deredundancy/0212/cross/OrthoFinder

mkdir -p step5_loci_ge80
for f in step4_loci_fa/*.fa; do
  n=$(grep -c '^>' "$f")
  [ "$n" -ge 54 ] && cp "$f" step5_loci_ge80/
done

ls step5_loci_ge80/*.fa | wc -l
