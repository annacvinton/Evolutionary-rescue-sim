#!/bin/bash
# ---------------------------------------------------------------------------
# Full design sweep. Resumable: re-running skips cells that already finished.
# Usage:   ./run_design.sh [n_parallel]     (default: all cores)
# ---------------------------------------------------------------------------
set -u
NPROC=${1:-$(getconf _NPROCESSORS_ONLN)}
BATCH=${BATCH:-}   # set BATCH=2 to add an independent batch of replicates
OUT=${OUT:-results_v2};  mkdir -p $OUT
JOBS=jobs_v2.txt

# ---- design ----------------------------------------------------------------
# All of these can be overridden from the environment, e.g.
#   PERTS="10 11" DRAWS=0 NREPS=5 ./run_design.sh
SLOPES=${SLOPES:-"0.8 1.0 1.2"}
DISPS=${DISPS:-"1.5 3 6"}
MUTS=${MUTS:-"0 0.75"}
PERTS=${PERTS:-"9 10 11 12 13"}
DRAWS=${DRAWS:-"0 1 2"}        # independent landscape realisations per condition
NREPS=${NREPS:-10}             # replicates per draw -> 30 per cell
TPERT=${TPERT:-250}            # perturbation time (equilibrates by ~t=60)
TMAX=${TMAX:-900}              # 650 post-shift, matching the original window
NFOUND=${NFOUND:-700}

# landscape conditions: SD 0 has no autocorrelation to speak of, so it is one
# condition, not three.  SD 1 and SD 3 each cross AC 0/2/4.  Seven in total.
LSCAPES=${LSCAPES:-"0:0 1:0 1:2 1:4 2:0 2:2 2:4"}

# ---- build the job list ----------------------------------------------------
: > $JOBS
for S in $SLOPES; do for D in $DISPS; do for M in $MUTS; do
for P in $PERTS; do for LS in $LSCAPES; do
  SD=${LS%%:*}; AC=${LS##*:}
  for R in $DRAWS; do
    TAG="S${S}_D${D}_M${M}_P${P}_sd${SD}_ac${AC}_r${R}${BATCH:+_b$BATCH}"
    [ -f "$OUT/${TAG}.done" ] && continue
    echo "$S $D $M $P $SD $AC $R $TAG" >> $JOBS
  done
done; done; done; done; done

TOTAL=$(wc -l < $JOBS)
RUNS=$((TOTAL*NREPS))
echo "$TOTAL cells queued = $RUNS runs on $NPROC cores"
echo "rough estimate: $((RUNS*27/10/60/NPROC)) hours"
echo "started $(date)"
[ -n "${DRYRUN:-}" ] && { echo "DRYRUN - nothing executed"; exit 0; }
[ "$TOTAL" -eq 0 ] && { echo "nothing to do"; exit 0; }

# ---- run -------------------------------------------------------------------
# One line per cell in $JOBS, eight fields each. runjob.sh receives them as
# positional arguments, so nothing long is ever constructed on the command line.
# (BSD xargs caps -I constructed arguments at 255 bytes, which the previous
# approach exceeded once the batch suffix was added.)
cat > runjob.sh <<'RUNNER'
#!/bin/sh
S=$1; D=$2; M=$3; P=$4; SD=$5; AC=$6; R=$7; TAG=$8
if [ "$SD" = "0" ]; then L="landscapes_v2/L_ac0_sd0_r${R}.txt"
else                     L="landscapes_v2/L_ac${AC}_sd${SD}_r${R}.txt"; fi
LANDSCAPE=$L OUTFILE=$OUT/${TAG}.csv SLOPE=$S DISP=$D MUTSD=$M PERT=$P \
  NREPS=$NREPS TMAX=$TMAX TPERT=$TPERT NFOUND=$NFOUND SEED=$$ \
  ./iiasa >/dev/null 2>&1 && touch $OUT/${TAG}.done
RUNNER
chmod +x runjob.sh
export NREPS TMAX TPERT NFOUND OUT

xargs -P "$NPROC" -L 1 ./runjob.sh < $JOBS
echo "complete $(date) -- $(ls $OUT/*.csv 2>/dev/null | wc -l) files in $OUT/"
