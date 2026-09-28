#!/bin/bash
# v3 CONTROL ARM: every treatment cell with NO perturbation, run to t=1500.
# Tests viability directly and gives the true unperturbed equilibrium per cell.
# 3 slopes x 3 disp x 2 mut x 7 landscapes = 126 cells x 3 draws x 3 reps = 1,134 runs.
set -u
NPROC=${1:-$(getconf _NPROCESSORS_ONLN)}
mkdir -p controls_v3
cat > runctrl.sh <<'RUNNER'
#!/bin/sh
S=$1; D=$2; M=$3; SD=$4; AC=$5; R=$6
TAG=C_S${S}_D${D}_M${M}_sd${SD}_ac${AC}_r${R}
[ -f controls_v3/${TAG}.done ] && exit 0
if [ "$SD" = "0" ]; then L=landscapes_v3/L_ac0_sd0_r${R}.txt; else L=landscapes_v3/L_ac${AC}_sd${SD}_r${R}.txt; fi
LANDSCAPE=$L OUTFILE=controls_v3/${TAG}.csv SLOPE=$S DISP=$D MUTSD=$M PERT=0 \
  NREPS=3 TMAX=1500 TPERT=99999 NFOUND=700 SEED=$$ ./iiasa_v3 >/dev/null 2>&1 \
  && touch controls_v3/${TAG}.done
RUNNER
chmod +x runctrl.sh
J=jobs_ctrl.txt; : > $J
for S in 0.8 1.0 1.2; do for D in 1.5 3 6; do for M in 0 0.75; do
  for SD in 0 1 2; do
    if [ "$SD" = "0" ]; then ACS="0"; else ACS="0 2 4"; fi
    for AC in $ACS; do for R in 0 1 2; do echo "$S $D $M $SD $AC $R" >> $J; done; done
  done
done; done; done
echo "$(wc -l < $J | tr -d ' ') cells x 3 reps, unperturbed to t=1500, on $NPROC cores   started $(date)"
xargs -P "$NPROC" -L 1 ./runctrl.sh < $J
echo "slope,disp,mutsd,patch_sd,ac,draw,rep,t,n" > controls_v3.csv
for f in controls_v3/C_*.csv; do
  b=$(basename "$f" .csv); b=${b#C_}
  p=$(echo "$b" | sed 's/^S//; s/_D/,/; s/_M/,/; s/_sd/,/; s/_ac/,/; s/_r/,/')
  awk -F, -v pre="$p" '($4==400||$4==800||$4==1200||$4==1500){print pre","$3","$4","$5}' "$f" >> controls_v3.csv
done
echo "complete $(date) -- $(($(wc -l < controls_v3.csv)-1)) rows in controls_v3.csv"
