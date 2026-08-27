#!/usr/bin/env bash
# MLST typing for the lineage-adjustment analysis (Reviewer 1, major concern 2).
#
# Batched download-then-type, so partial progress survives an interruption and so
# peak disk stays bounded: each batch's archive is deleted once extracted.
#
# The scheme is FORCED per organism rather than auto-detected. A. baumannii has two
# schemes in PubMLST (Oxford and Pasteur) and letting mlst choose per genome would
# mix them, producing STs that are not comparable to one another. Pasteur
# (abaumannii_2) is used because it is the scheme with stable locus definitions.
#
# Every accession ends with exactly one status row. A genome that produced no FASTA,
# or no ST line, is recorded as a failure - never skipped silently.

set -u
cd "$HOME/mlst_run"
E="$HOME/micromamba/envs/mlst"
export PATH="$E/bin:$PATH"

BATCH=100
IN=mlst_accessions.tsv
OUT=mlst_calls.tsv
STATUS=mlst_status.tsv
LOG=mlst_run.log
THREADS=100

mkdir -p genomes work
[ -f "$OUT" ]    || printf 'assembly_version\torganism\tscheme\tst\tprofile\n' > "$OUT"
[ -f "$STATUS" ] || printf 'assembly_version\tstage\tstatus\tdetail\n' > "$STATUS"

done_acc() { cut -f1 "$STATUS" | tail -n +2 | sort -u; }

scheme_for() {
  case "$1" in
    "Acinetobacter baumannii") echo abaumannii_2 ;;
    Klebsiella*)               echo klebsiella ;;
    *)                         echo "" ;;
  esac
}

echo "=== start $(date -u +%FT%TZ) ===" >> "$LOG"
tail -n +2 "$IN" > work/all.tsv
done_acc > work/done.txt
# An empty first file makes NR==FNR true for the whole of the second file, which
# would drop every accession and print "to process: 0" as though the run were
# already finished. That failure looks exactly like success, so it is guarded.
if [ -s work/done.txt ]; then
  awk -F'\t' 'NR==FNR{d[$1];next} !($1 in d)' work/done.txt work/all.tsv > work/todo.tsv
else
  cp work/all.tsv work/todo.tsv
fi
TOTAL=$(wc -l < work/todo.tsv)
echo "to process: $TOTAL" | tee -a "$LOG"

split -l "$BATCH" -d -a 4 work/todo.tsv work/b_
for B in work/b_*; do
  n=$(wc -l < "$B")
  echo "--- batch $(basename "$B") ($n) $(date -u +%T)" >> "$LOG"
  cut -f1 "$B" > work/acc.txt
  rm -rf work/dl work/pkg.zip
  if ! datasets download genome accession --inputfile work/acc.txt \
        --include genome --no-progressbar --filename work/pkg.zip >> "$LOG" 2>&1; then
    while read -r a o; do
      printf '%s\tdownload\tFAIL\tdatasets returned non-zero\n' "$a" >> "$STATUS"
    done < "$B"
    continue
  fi
  unzip -qq -o work/pkg.zip -d work/dl >> "$LOG" 2>&1
  rm -f work/pkg.zip

  while IFS=$'\t' read -r acc org; do
    f=$(find work/dl -type f -name "${acc}*_genomic.fna" 2>/dev/null | head -1)
    if [ -z "$f" ]; then
      f=$(find work/dl/ncbi_dataset/data/"$acc" -type f -name "*.fna" 2>/dev/null | head -1)
    fi
    if [ -z "$f" ] || [ ! -s "$f" ]; then
      printf '%s\tdownload\tFAIL\tno FASTA in the package\n' "$acc" >> "$STATUS"
      continue
    fi
    cp "$f" "genomes/${acc}.fna"
    printf '%s\tdownload\tOK\t%s bytes\n' "$acc" "$(stat -c%s "genomes/${acc}.fna")" >> "$STATUS"
  done < "$B"
  rm -rf work/dl

  # ---- type this batch in parallel -----------------------------------------
  : > work/batch_calls.tsv
  while IFS=$'\t' read -r acc org; do
    [ -s "genomes/${acc}.fna" ] || continue
    s=$(scheme_for "$org")
    [ -n "$s" ] || { printf '%s\tmlst\tFAIL\tno scheme registered for %s\n' "$acc" "$org" >> "$STATUS"; continue; }
    printf '%s\t%s\t%s\n' "$acc" "$org" "$s"
  done < "$B" > work/batch_todo.tsv

  cat work/batch_todo.tsv | xargs -P "$THREADS" -I{} -d'\n' bash -c '
    IFS="'"$(printf '\t')"'" read -r acc org sch <<< "{}"
    line=$(mlst --quiet --nopath --scheme "$sch" "genomes/${acc}.fna" 2>/dev/null | head -1)
    if [ -z "$line" ]; then
      printf "%s\tmlst\tFAIL\tno output line\n" "$acc" >> mlst_status.tsv
    else
      st=$(printf "%s" "$line" | cut -f3)
      prof=$(printf "%s" "$line" | cut -f4- | tr "\t" ",")
      printf "%s\t%s\t%s\t%s\t%s\n" "$acc" "$org" "$sch" "$st" "$prof" >> mlst_calls.tsv
      printf "%s\tmlst\tOK\tST %s\n" "$acc" "$st" >> mlst_status.tsv
    fi
  '
  rm -f genomes/*.fna
  echo "  batch done: calls=$(( $(wc -l < "$OUT") - 1 ))" >> "$LOG"
done

echo "=== finished $(date -u +%FT%TZ) ===" >> "$LOG"
echo "TOTAL_CALLS=$(( $(wc -l < "$OUT") - 1 ))" >> "$LOG"
echo "MLST_RUN_COMPLETE" >> "$LOG"
