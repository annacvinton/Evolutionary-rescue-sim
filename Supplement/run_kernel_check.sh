#!/bin/bash
# Kernel sensitivity check: exponential (original) vs Gaussian competition kernel.
# Same treatment cell as the trajectory figure -- slope 1.0, dispersal 3, mutation on,
# severity 11 -- on the flat landscape and patch SD 2 at AC 0/2/4. 3 draws x 4 reps
# per condition per kernel = 96 runs. Full v2 timing (perturb t=250, run to t=900).
set -u
NPROC=${1:-$(getconf _NPROCESSORS_ONLN)}
mkdir -p results_kernel
cat > runkernel.sh <<'RUNNER'
#!/bin/sh
K=$1; SD=$2; AC=$3; R=$4; TAG=k${K}_sd${SD}_ac${AC}_r${R}
[ -f results_kernel/${TAG}.done ] && exit 0
if [ "$SD" = "0" ]; then L=landscapes_v2/L_ac0_sd0_r${R}.txt; else L=landscapes_v2/L_ac${AC}_sd${SD}_r${R}.txt; fi
KERNEL=$K LANDSCAPE=$L OUTFILE=results_kernel/${TAG}.csv SLOPE=1.0 DISP=3 MUTSD=0.75 PERT=11 \
  NREPS=4 TMAX=900 TPERT=250 NFOUND=700 SEED=$$ ./iiasa_kernel >/dev/null 2>&1 \
  && touch results_kernel/${TAG}.done
RUNNER
chmod +x runkernel.sh
: > jobs_kernel.txt
for K in 0 1; do for cond in "0 0" "2 0" "2 2" "2 4"; do set -- $cond; for R in 0 1 2; do echo "$K $1 $2 $R" >> jobs_kernel.txt; done; done; done
echo "$(wc -l < jobs_kernel.txt) cells on $NPROC cores"
xargs -P "$NPROC" -L 1 ./runkernel.sh < jobs_kernel.txt
echo "kernel,patch_sd,ac,draw,rep,t,n" > kernel_check.csv
for f in results_kernel/k*.csv; do
  b=$(basename "$f" .csv); K=${b:1:1}; rest=${b#*_sd}; SD=${rest%%_*}; rest=${rest#*_ac}; AC=${rest%%_*}; R=${rest##*_r}
  awk -F, -v k=$K -v sd=$SD -v ac=$AC -v r=$R '{print k","sd","ac","r","$3","$4","$5}' "$f" >> kernel_check.csv
done
echo "complete $(date) -- $(($(wc -l < kernel_check.csv)-1)) rows in kernel_check.csv"
