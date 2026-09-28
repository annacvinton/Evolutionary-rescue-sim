#!/bin/bash
# v3 PILOT. Perturbation at t=400 (run to 1050). Gaussian kernel, founders matched to
# local optimum, RNG tied to SEED, y-periodic competition, maladaptation includes p(t).
# Adults buffered by default (ADULTS_SHIFT=0); a few ADULTS_SHIFT=1 runs at low severity
# are included for the record. Slope 0 (homogeneous) on the flat landscape only.
set -u
NPROC=${1:-$(getconf _NPROCESSORS_ONLN)}
mkdir -p results_pilot_v3
cat > runpilot.sh <<'RUNNER'
#!/bin/sh
A=$1; S=$2; D=$3; M=$4; P=$5; SD=$6; AC=$7; R=$8
TAG=A${A}_S${S}_D${D}_M${M}_P${P}_sd${SD}_ac${AC}_r${R}
[ -f results_pilot_v3/${TAG}.done ] && exit 0
if [ "$SD" = "0" ]; then L=landscapes_v2/L_ac0_sd0_r${R}.txt; else L=landscapes_v2/L_ac${AC}_sd${SD}_r${R}.txt; fi
ADULTS_SHIFT=$A LANDSCAPE=$L OUTFILE=results_pilot_v3/${TAG}.csv SLOPE=$S DISP=$D MUTSD=$M PERT=$P \
  NREPS=4 TMAX=1050 TPERT=400 NFOUND=700 SEED=$$ ./iiasa_v3 >/dev/null 2>&1 \
  && touch results_pilot_v3/${TAG}.done
RUNNER
chmod +x runpilot.sh
J=jobs_pilot_v3.txt; : > $J
L3="0 0
2 0
2 4"
# core: slope x severity, disp 3, mut on
for S in 0.8 1.0 1.2; do for P in 9 11 13; do echo "$L3" | while read SD AC; do echo "0 $S 3 0.75 $P $SD $AC 0" >> $J; done; done; done
# dispersal extremes
for D in 1.5 6; do for P in 9 11 13; do echo "$L3" | while read SD AC; do echo "0 1.0 $D 0.75 $P $SD $AC 0" >> $J; done; done; done
# mutation off
for P in 9 11 13; do echo "$L3" | while read SD AC; do echo "0 1.0 3 0 $P $SD $AC 0" >> $J; done; done
# homogeneous (slope 0), flat landscape only -- EXPENSIVE (N~3600)
for M in 0 0.75; do for P in 9 13; do echo "0 0 3 $M $P 0 0 0" >> $J; done; done
# adults feel the shift, low severities, flat -- for the record
for P in 3 5 7 9; do echo "1 1.0 3 0.75 $P 0 0 0" >> $J; done
echo "$(wc -l < $J) cells x 4 reps on $NPROC cores   started $(date)"
xargs -P "$NPROC" -L 1 ./runpilot.sh < $J
echo "adults_shift,slope,disp,mutsd,pert,patch_sd,ac,draw,rep,t,n,u_mean,u_sd,mal_mean,x_sd" > pilot_v3.csv
for f in results_pilot_v3/A*.csv; do
  b=$(basename "$f" .csv)
  p=$(echo "$b" | sed 's/^A//; s/_S/,/; s/_D/,/; s/_M/,/; s/_P/,/; s/_sd/,/; s/_ac/,/; s/_r/,/')
  awk -F, -v pre="$p" '{print pre","$3","$4","$5","$6","$7","$10","$19}' "$f" >> pilot_v3.csv
done
echo "complete $(date) -- $(($(wc -l < pilot_v3.csv)-1)) rows in pilot_v3.csv"
