#!/usr/bin/env bash
# NM-V1 independent parser verification. Reads RAW tool outputs only.
# Shares no code with the collection script or the Python builder: pure awk, different logic.
set -uo pipefail
# Both directories are required and deliberately have no default: the
# original values named a specific machine and user account.
R="${NMV1_RUN_DIR:?set NMV1_RUN_DIR to the directory holding the raw tool outputs}"
C="${NMV1_DERIVED_DIR:?set NMV1_DERIVED_DIR to the collected nmv1 derived directory}"
FAIL=0
chk () { # name expected actual
  if [ "$2" = "$3" ]; then printf "  %-42s expected %-6s got %-6s MATCH\n" "$1" "$2" "$3";
  else printf "  %-42s expected %-6s got %-6s *** MISMATCH ***\n" "$1" "$2" "$3"; FAIL=$((FAIL+1)); fi
}

echo "=== A. ISEScan, recounted from raw .fna.tsv files ==="
# total IS elements: every non-header non-empty row across all raw tsv
TOT=$(find "$R/isescan_out/fasta" -name '*.fna.tsv' -exec awk 'FNR>1 && NF>0' {} + | wc -l)
# complete / partial from column 22 of the RAW file
COMP=$(find "$R/isescan_out/fasta" -name '*.fna.tsv' -exec awk -F'\t' 'FNR>1 && $22=="c"' {} + | wc -l)
PART=$(find "$R/isescan_out/fasta" -name '*.fna.tsv' -exec awk -F'\t' 'FNR>1 && $22=="p"' {} + | wc -l)
chk "total IS elements"            1648 "$TOT"
chk "complete IS elements"         1215 "$COMP"
chk "partial IS elements"          433  "$PART"
chk "complete + partial == total"  "$TOT" "$(( COMP + PART ))"

POSB=0; COMPB=0; PARTONLY=0
for f in "$R"/isescan_out/fasta/*.fna.tsv; do
  [ -e "$f" ] || continue
  n=$(awk 'FNR>1 && NF>0' "$f" | wc -l)
  [ "$n" -gt 0 ] && POSB=$((POSB+1))
  c=$(awk -F'\t' 'FNR>1 && $22=="c"' "$f" | wc -l)
  p=$(awk -F'\t' 'FNR>1 && $22=="p"' "$f" | wc -l)
  [ "$c" -gt 0 ] && COMPB=$((COMPB+1))
  [ "$c" -eq 0 ] && [ "$p" -gt 0 ] && PARTONLY=$((PARTONLY+1))
done
chk "IS-positive blocks"           602 "$POSB"
chk "blocks with >=1 complete IS"  491 "$COMPB"
chk "partial-only blocks"          111 "$PARTONLY"

echo
echo "=== B. IntegronFinder, recounted from raw .summary files ==="
NSUM=$(find "$R/if_out" -name '*.summary' | wc -l)
chk "summary files present"        1282 "$NSUM"
chk "missing IntegronFinder output" 1 "$(( 1283 - NSUM ))"
IPOS=0; ICOMP=0
for f in $(find "$R/if_out" -name '*.summary'); do
  read ca ci i0 <<< "$(grep -v '^#' "$f" | awk -F'\t' 'NR>1 && NF>=4 {a+=$2; c+=$3; z+=$4} END{print a+0,c+0,z+0}')"
  [ $(( ca + ci + i0 )) -gt 0 ] && IPOS=$((IPOS+1))
  [ "$ci" -gt 0 ] && ICOMP=$((ICOMP+1))
done
chk "integron-positive blocks"     196 "$IPOS"
chk "complete-integron blocks"     163 "$ICOMP"

echo
echo "=== C. cross-check the COLLECTED tables against this recount ==="
chk "collected isescan_hits rows"  "$TOT" "$(( $(wc -l < "$C/isescan_hits.tsv") - 1 ))"
chk "collected per-block rows"     1283 "$(( $(wc -l < "$C/isescan_per_block.tsv") - 1 ))"
chk "collected IF per-block rows"  1283 "$(( $(wc -l < "$C/integronfinder_per_block.tsv") - 1 ))"
chk "collected hits type=c"        "$COMP" "$(awk -F'\t' 'NR>1 && $23=="c"' "$C/isescan_hits.tsv" | wc -l)"
chk "collected hits type=p"        "$PART" "$(awk -F'\t' 'NR>1 && $23=="p"' "$C/isescan_hits.tsv" | wc -l)"
chk "collected IS-positive blocks" "$POSB" "$(awk -F'\t' 'NR>1 && $5=="yes"' "$C/isescan_per_block.tsv" | wc -l)"
chk "collected complete-IS blocks" "$COMPB" "$(awk -F'\t' 'NR>1 && $3>0' "$C/isescan_per_block.tsv" | wc -l)"
chk "collected partial-only"       "$PARTONLY" "$(awk -F'\t' 'NR>1 && $3==0 && $4>0' "$C/isescan_per_block.tsv" | wc -l)"
chk "collected IF-positive blocks" "$IPOS" "$(awk -F'\t' 'NR>1 && $5=="yes"' "$C/integronfinder_per_block.tsv" | wc -l)"
chk "collected IF complete blocks" "$ICOMP" "$(awk -F'\t' 'NR>1 && $3>0' "$C/integronfinder_per_block.tsv" | wc -l)"
chk "collected IF missing_output"  1 "$(awk -F'\t' 'NR>1 && $6!="ok"' "$C/integronfinder_per_block.tsv" | wc -l)"

echo
echo "=== D. per-block field agreement, all 1283 blocks ==="
BAD=0
while IFS=$'\t' read -r b n c p pos; do
  [ "$b" = "block_id" ] && continue
  t="$R/isescan_out/fasta/$b.fna.tsv"
  if [ -f "$t" ]; then
    rn=$(awk 'FNR>1 && NF>0' "$t" | wc -l)
    rc=$(awk -F'\t' 'FNR>1 && $22=="c"' "$t" | wc -l)
    rp=$(awk -F'\t' 'FNR>1 && $22=="p"' "$t" | wc -l)
  else rn=0; rc=0; rp=0; fi
  if [ "$rn" != "$n" ] || [ "$rc" != "$c" ] || [ "$rp" != "$p" ]; then
    BAD=$((BAD+1)); [ "$BAD" -le 3 ] && echo "     MISMATCH $b: raw($rn,$rc,$rp) vs collected($n,$c,$p)"
  fi
done < "$C/isescan_per_block.tsv"
chk "per-block IS fields disagreeing" 0 "$BAD"

echo
if [ "$FAIL" -eq 0 ]; then echo "VERIFIER VERDICT: PASS - zero disagreements"; else echo "VERIFIER VERDICT: *** FAIL - $FAIL disagreements ***"; fi
exit $FAIL
