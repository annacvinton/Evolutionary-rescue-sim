#!/bin/bash
# Kernel sensitivity, extended: one-factor-at-a-time star around the centre cell
# (slope 1.0, dispersal 3, mutation 0.75, severity 11), on flat + patch SD 2 at
# AC 0/2/4, exponential vs Gaussian kernel. 3 draws x 3 reps per cell per kernel.
# 8 treatment cells x 4 landscapes x 3 draws x 2 kernels = 192 jobs, 576 runs.
set -u
NPROC=${1:-$(getconf _NPROCESSORS_ONLN)}
mkdir -p results_kernel2
cat > runkernel2.sh <<'RUNNER'
#!/bin/sh
K=$1; S=$2; D=$3; M=$4; P=$5; SD=$6; AC=$7; R=$8
TAG=k${K}_S${S}_D${D}_M${M}_P${P}_sd${SD}_ac${AC}_r${R}
[ -f results_kernel2/${TAG}.done ] && exit 0
if [ "$SD" = "0" ]; then L=landscapes_v2/L_ac0_sd0_r${R}.txt; else L=landscapes_v2/L_ac${AC}_sd${SD}_r${R}.txt; fi
KERNEL=$K LANDSCAPE=$L OUTFILE=results_kernel2/${TAG}.csv SLOPE=$S DISP=$D MUTSD=$M PERT=$P \
  NREPS=3 TMAX=900 TPERT=250 NFOUND=700 SEED=$$ ./iiasa_kernel >/dev/null 2>&1 \
  && touch results_kernel2/${TAG}.done
RUNNER
chmod +x runkernel2.sh
: > jobs_kernel2.txt
# centre + one-factor-at-a-time extremes
CELLS="1.0 3 0.75 11
0.8 3 0.75 11
1.2 3 0.75 11
1.0 1.5 0.75 11
1.0 6 0.75 11
1.0 3 0 11
1.0 3 0.75 9
1.0 3 0.75 13"
echo "$CELLS" | while read S D M P; do
  for K in 0 1; do for cond in "0 0" "2 0" "2 2" "2 4"; do set -- $cond; for R in 0 1 2; do
    echo "$K $S $D $M $P $1 $2 $R" >> jobs_kernel2.txt
  done; done; done
done
echo "$(wc -l < jobs_kernel2.txt) jobs on $NPROC cores"
xargs -P "$NPROC" -L 1 ./runkernel2.sh < jobs_kernel2.txt
echo "kernel,slope,disp,mutsd,pert,patch_sd,ac,draw,rep,t,n" > kernel_check2.csv
for f in results_kernel2/k*.csv; do
  b=$(basename "$f" .csv)
  p=$(echo "$b" | sed 's/^k//; s/_S/,/; s/_D/,/; s/_M/,/; s/_P/,/; s/_sd/,/; s/_ac/,/; s/_r/,/')
  awk -F, -v pre="$p" '{print pre","$3","$4","$5}' "$f" >> kernel_check2.csv
done
echo "complete $(date) -- $(($(wc -l < kernel_check2.csv)-1)) rows in kernel_check2.csv"
