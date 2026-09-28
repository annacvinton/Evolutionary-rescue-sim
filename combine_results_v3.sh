#!/bin/bash
# Concatenate v3 per-cell CSVs into all_results_v3.csv with treatment columns.
# v3 filenames: S{slope}_D{disp}_M{mut}_P{pert}_sd{SD}_ac{AC}_r{draw}.csv  (no batch; batch=1)
set -u
echo "slope,disp,mutsd,pert_treat,patch_sd,ac,draw,batch,pert_value,pert_name,rep,t,n,u_mean,u_sd,u_skew,u_kurt,mal_mean,mal_sd,mal_skew,mal_kurt,nn_mean,nn_sd,nn_skew,nn_kurt,x_mean,x_sd,x_skew,x_kurt" > all_results_v3.csv
for fp in results_v3/S*.csv; do
  b=$(basename "$fp" .csv)
  p=$(echo "$b" | sed 's/^S//; s/_D/,/; s/_M/,/; s/_P/,/; s/_sd/,/; s/_ac/,/; s/_r/,/; s/$/,1/')
  awk -v pre="$p" -F, 'NF==21{print pre","$0}' "$fp" >> all_results_v3.csv
done
echo "$(ls results_v3/S*.csv | wc -l | tr -d ' ') files -> all_results_v3.csv ($(($(wc -l < all_results_v3.csv)-1)) rows)"
