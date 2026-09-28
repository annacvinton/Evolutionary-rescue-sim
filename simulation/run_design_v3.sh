#!/bin/bash
# ============================================================================
# v3 FULL SWEEP
#   630 treatment combinations x 5 landscape draws x 10 runs = 31,500 runs
#   + homogeneous block: slope 0, flat, disp 3, mut {0,0.75}, pert 9-13,
#     10 runs each = 100 runs
#   Perturbation t=400, run to t=1050. Buffered adults (ADULTS_SHIFT=0).
#   Resumable: completed cells carry a .done marker and are skipped.
#   Usage:  ./run_design_v3.sh [ncores]
# ============================================================================
set -u
NPROC=${1:-$(getconf _NPROCESSORS_ONLN)}
mkdir -p results_v3
cat > runone_v3.sh <<'RUNNER'
#!/bin/sh
S=$1; D=$2; M=$3; P=$4; SD=$5; AC=$6; R=$7
TAG=S${S}_D${D}_M${M}_P${P}_sd${SD}_ac${AC}_r${R}
[ -f results_v3/${TAG}.done ] && exit 0
if [ "$SD" = "0" ]; then L=landscapes_v3/L_ac0_sd0_r${R}.txt; else L=landscapes_v3/L_ac${AC}_sd${SD}_r${R}.txt; fi
LANDSCAPE=$L OUTFILE=results_v3/${TAG}.csv SLOPE=$S DISP=$D MUTSD=$M PERT=$P \
  NREPS=10 TMAX=1050 TPERT=400 NFOUND=700 SEED=$$ ./iiasa_v3 >/dev/null 2>&1 \
  && touch results_v3/${TAG}.done
RUNNER
chmod +x runone_v3.sh
J=jobs_v3.txt; : > $J
for S in 0.8 1.0 1.2; do
 for D in 1.5 3 6; do
  for M in 0 0.75; do
   for P in 9 10 11 12 13; do
    for SD in 0 1 2; do
     if [ "$SD" = "0" ]; then ACS="0"; else ACS="0 2 4"; fi
     for AC in $ACS; do
      for R in 0 1 2 3 4; do echo "$S $D $M $P $SD $AC $R" >> $J; done
     done
    done
   done
  done
 done
done
# homogeneous block: slope 0, flat only, disp 3
for M in 0 0.75; do for P in 9 10 11 12 13; do echo "0 3 $M $P 0 0 0" >> $J; done; done
N=$(wc -l < $J | tr -d ' ')
echo "$N cells queued = $((N*10)) runs on $NPROC cores"
echo "rough estimate: 45-50 hours; slope-0 cells are the slow tail"
echo "started $(date)"
xargs -P "$NPROC" -L 1 ./runone_v3.sh < $J
echo "complete $(date) -- $(ls results_v3/*.done 2>/dev/null | wc -l) cells done"
