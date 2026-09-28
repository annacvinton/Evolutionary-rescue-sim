#!/bin/bash
# Slope 0.6 block for v3: completes the balanced factorial 0.6/0.8/1.0.
#   Sweep: slope 0.6 x 3 disp x 2 mut x 5 pert x 7 landscapes x 5 draws x 10 runs
#          = 1,050 cells = 10,500 runs, into results_v3/ (same naming; combine
#          and reduce scripts pick them up automatically)
#   Controls: slope 0.6 x 3 disp x 2 mut x 7 landscapes x 3 draws, unperturbed
#          to t=1500, 3 reps = 126 x ... -> controls_v3/ (same naming)
# Resumable. Usage: ./run_slope06_v3.sh [ncores]
set -u
NPROC=${1:-$(getconf _NPROCESSORS_ONLN)}
mkdir -p results_v3 controls_v3
cat > runone06.sh <<'RUNNER'
#!/bin/sh
MODE=$1; S=0.6; D=$2; M=$3; P=$4; SD=$5; AC=$6; R=$7
if [ "$SD" = "0" ]; then L=landscapes_v3/L_ac0_sd0_r${R}.txt; else L=landscapes_v3/L_ac${AC}_sd${SD}_r${R}.txt; fi
if [ "$MODE" = "sweep" ]; then
  TAG=S${S}_D${D}_M${M}_P${P}_sd${SD}_ac${AC}_r${R}
  [ -f results_v3/${TAG}.done ] && exit 0
  LANDSCAPE=$L OUTFILE=results_v3/${TAG}.csv SLOPE=$S DISP=$D MUTSD=$M PERT=$P \
    NREPS=10 TMAX=1050 TPERT=400 NFOUND=700 SEED=$$ ./iiasa_v3 >/dev/null 2>&1 \
    && touch results_v3/${TAG}.done
else
  TAG=C_S${S}_D${D}_M${M}_sd${SD}_ac${AC}_r${R}
  [ -f controls_v3/${TAG}.done ] && exit 0
  LANDSCAPE=$L OUTFILE=controls_v3/${TAG}.csv SLOPE=$S DISP=$D MUTSD=$M PERT=0 \
    NREPS=3 TMAX=1500 TPERT=99999 NFOUND=700 SEED=$$ ./iiasa_v3 >/dev/null 2>&1 \
    && touch controls_v3/${TAG}.done
fi
RUNNER
chmod +x runone06.sh
J=jobs_s06.txt; : > $J
for D in 1.5 3 6; do for M in 0 0.75; do
  for SD in 0 1 2; do
    if [ "$SD" = "0" ]; then ACS="0"; else ACS="0 2 4"; fi
    for AC in $ACS; do
      for R in 0 1 2; do echo "ctrl $D $M 0 $SD $AC $R" >> $J; done
      for P in 9 10 11 12 13; do for R in 0 1 2 3 4; do echo "sweep $D $M $P $SD $AC $R" >> $J; done; done
    done
  done
done; done
echo "$(grep -c ^ctrl $J) control cells + $(grep -c ^sweep $J) sweep cells on $NPROC cores   started $(date)"
xargs -P "$NPROC" -L 1 ./runone06.sh < $J
echo "complete $(date)"
echo "next: rerun ./combine_results_v3.sh and python3 reduce_sweep_v3.py (they include slope 0.6 automatically),"
echo "and re-collect controls: bash -c 'tail -n +1 /dev/null'; see chat for the control re-collect line"
