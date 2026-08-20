#!/bin/bash
# Supplementary viability map: unperturbed runs across patch SD x autocorrelation.
# 30 cells x 3 draws x 4 reps = 360 runs, no perturbation, t=800.
# Resumable; writes boundary_map.csv incrementally.
set -u
NPROC=${1:-$(getconf _NPROCESSORS_ONLN)}
mkdir -p results_boundary
[ -f boundary_map.csv ] || echo "sd,ac,draw,rep,final_t,final_n" > boundary_map.csv

cat > runboundary.sh <<'RUNNER'
#!/bin/sh
SD=$1; AC=$2; R=$3; TAG=bnd_${SD}_${AC}_${R}
[ -f results_boundary/${TAG}.done ] && exit 0
LANDSCAPE=landscapes_boundary/B_ac${AC}_sd${SD}_r${R}.txt \
  OUTFILE=results_boundary/${TAG}.csv SLOPE=1.0 DISP=3 MUTSD=0.75 PERT=0 \
  NREPS=4 TMAX=800 TPERT=99999 NFOUND=700 SEED=$$ ./iiasa >/dev/null 2>&1 \
  && touch results_boundary/${TAG}.done
RUNNER
chmod +x runboundary.sh

JOBS=jobs_boundary.txt; : > $JOBS
for SD in 0.5 1.0 1.5 2.0 2.5 3.0; do
  for AC in 0 1 2 3 4; do
    for R in 0 1 2; do echo "$SD $AC $R" >> $JOBS; done
  done
done
echo "$(wc -l < $JOBS) cells on $NPROC cores"
xargs -P "$NPROC" -L 1 ./runboundary.sh < $JOBS

# collect
echo "sd,ac,draw,rep,final_t,final_n" > boundary_map.csv
for f in results_boundary/bnd_*.csv; do
  b=$(basename "$f" .csv); b=${b#bnd_}
  SD=${b%%_*}; rest=${b#*_}; AC=${rest%%_*}; R=${rest##*_}
  awk -F, -v sd=$SD -v ac=$AC -v r=$R '{k=$3;t[k]=$4;n[k]=$5}
    END{for(x in t) print sd","ac","r","x","t[x]","n[x]}' "$f" >> boundary_map.csv
done
echo "complete $(date) -- $(($(wc -l < boundary_map.csv)-1)) runs in boundary_map.csv"
