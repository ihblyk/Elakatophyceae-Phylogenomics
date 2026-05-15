cd /mnt/sdb/geminella/cross/OrthoFinder/Results_May11
mkdir -p step6_trim_renamed

for f in step6_trim/*.trim.fa; do
  b=$(basename "$f")
  awk '
    /^>/ {
      sub(/^>/,"",$0)
      split($0,a,"|")
      print ">" a[1]
      next
    }
    {print}
  ' "$f" > "step6_trim_renamed/$b"
done
grep -R "^>.*|" step6_trim_renamed | head
