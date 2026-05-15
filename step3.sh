cd /mnt/sdb/transcriptome/elakatophyceae_pep/pp/deredundancy/0212/cross
mkdir -p hmmsearch_tbl hmmsearch_domtbl hmmsearch_log

QUERY_DIR=/mnt/sdb/transcriptome/elakatophyceae_pep/pp/deredundancy/0212/cross/OrthoFinder/all

JOBS=100
find "$QUERY_DIR" -name "*.pep" | xargs -I{} -P ${JOBS} bash -c '
q="$1"
base=$(basename "$q" .pep)
hmmsearch --cpu 2 -E 1e-6 --domE 1e-6\
  --tblout "hmmsearch_tbl/${base}.tbl" \
  --domtblout "hmmsearch_domtbl/${base}.domtbl" \
  SOG468.hmmdb "$q" > "hmmsearch_log/${base}.log" 2>&1
' _ {}
