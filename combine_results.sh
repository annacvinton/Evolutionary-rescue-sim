#!/bin/bash
# =============================================================================
# Concatenate the per-cell simulation output into one table.
#
#   ./combine_results.sh          -> writes all_results_combined.csv
#
# Each file in results/ is one cell of the design, named
#   S{slope}_D{disp}_M{mutsd}_P{severity}_sd{patchSD}_ac{AC}_r{draw}[_b{batch}]
# and contains one row per recorded timestep with no treatment columns. This
# parses the treatment levels back out of the filename and prepends them, so
# the result is self-describing.
#
# Files without a _b suffix are batch 1; _b2 is the second batch of replicates.
# Batch has to be carried through because the two batches reuse replicate
# numbers, so (rep) alone does not identify a run.
# =============================================================================
set -u
OUT=${1:-all_results_combined.csv}

echo "slope,disp,mutsd,pert_treat,patch_sd,ac,draw,batch,pert_value,pert_name,rep,t,n,u_mean,u_sd,u_skew,u_kurt,mal_mean,mal_sd,mal_skew,mal_kurt,nn_mean,nn_sd,nn_skew,nn_kurt,x_mean,x_sd,x_skew,x_kurt" > "$OUT"

n=0
for f in results/*.csv; do
  b=$(basename "$f" .csv)
  case "$b" in
    *_b2) batch=2; b=${b%_b2} ;;
    *)    batch=1 ;;
  esac
  # S0.8_D1.5_M0_P9_sd3_ac4_r2  ->  0.8,1.5,0,9,3,4,2
  p=$(echo "$b" | sed "s/^S//; s/_D/,/; s/_M/,/; s/_P/,/; s/_sd/,/; s/_ac/,/; s/_r/,/; s/\$/,$batch/")
  awk -v pre="$p" '{print pre","$0}' "$f" >> "$OUT"
  n=$((n+1))
done

echo "combined $n files into $OUT"
echo "  header fields: $(head -1 "$OUT" | awk -F, '{print NF}')  (expect 29)"
echo "  data fields:   $(awk -F, 'NR==2{print NF; exit}' "$OUT")  (expect 29)"
echo "  rows:          $(( $(wc -l < "$OUT") - 1 ))"
