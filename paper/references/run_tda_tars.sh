#!/usr/bin/env bash
# Triplet-window convergence check for the paper on tars: the pipeline TDA window against a 20 eV window and the
# full TDA (run_tda_full.py). Results are pushed to the branch results-tars-tda as Parfait after every finished job
# (never force-pushed).
# Usage (inside a worktree of the paper branch, psscreen conda env active):
#   nohup bash paper/references/run_tda_tars.sh > paper/references/logs/tda.log 2>&1 &
# Jobs run JOBS at a time (default 3), each with THREADS OpenMP threads (default 16).
set -uo pipefail
export PYTHONUNBUFFERED=1
cd "$(dirname "$0")/../.."
THREADS="${THREADS:-16}"
JOBS="${JOBS:-3}"
MEM_MB="${PSS_MAX_MEMORY_MB:-60000}"
LOGDIR=paper/references/logs
mkdir -p "$LOGDIR" paper/references/results paper/quantum/results
python -c "import ase" 2>/dev/null || pip install -q ase

# Optional: pass job names as arguments to run only those (e.g. "tda_full_azaBDP_Br2").
JOBLIST=(
  "paper/references/run_tda_full.py azaBDP_Br2"
  "paper/references/run_tda_full.py I2-BDP-0"
  "paper/references/run_tda_full.py I2-BDP"
)

push() {
  git add paper/references/results/*.tda_full.json 2>/dev/null
  if git diff --cached --quiet; then return 0; fi
  git -c user.name="Parfait" -c user.email="26543855+pifparfait@users.noreply.github.com" \
    commit -q -m "Add full-TDA triplet check: $1 ($(date -u +%Y-%m-%d))"
  if ! git push -q origin "HEAD:refs/heads/results-tars-tda"; then
    alt="results-tars-tda-$(date -u +%Y%m%d-%H%M%S)"
    echo "Push rejected; not forcing. Publishing to $alt instead." >&2
    git push -q origin "HEAD:refs/heads/$alt"
  fi
}

run_one() {
  local tag; tag=$(echo "$1" | awk '{print $1"_"$2}' | sed 's|.*/||; s|\.py||')
  OMP_NUM_THREADS=$THREADS PSS_MAX_MEMORY_MB=$MEM_MB python $1 > "$LOGDIR/tda_${tag}.log" 2>&1
  echo "finished $1 (exit $?)"
  ( flock 9; push "$tag" ) 9>"/tmp/psscreen_tda_$USER.lock"   # parallel jobs push one at a time
}

for job in "${JOBLIST[@]}"; do
  tag=$(echo "$job" | awk '{print $1"_"$2}' | sed 's|.*/||; s|\.py||; s|^run_||')
  if [ "$#" -gt 0 ] && [[ " $* " != *" $tag "* ]]; then continue; fi
  while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do wait -n; done
  run_one "$job" &
  echo "started $job"
done
wait
push "all"
echo "all TDA jobs finished"
