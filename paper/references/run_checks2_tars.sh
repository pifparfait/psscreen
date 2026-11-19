#!/usr/bin/env bash
# Second batch of checks for the paper on tars: all-electron SOC benchmark, vertical-to-adiabatic corrections and the
# THC factorization error. Results are pushed to the branch results-tars-checks2 as Parfait after every finished job
# (never force-pushed).
# Usage (inside a worktree of the paper branch, psscreen conda env active):
#   nohup bash paper/references/run_checks2_tars.sh > paper/references/logs/checks2.log 2>&1 &
# Jobs run JOBS at a time (default 3), each with THREADS OpenMP threads (default 16).
set -uo pipefail
cd "$(dirname "$0")/../.."
THREADS="${THREADS:-16}"
JOBS="${JOBS:-3}"
MEM_MB="${PSS_MAX_MEMORY_MB:-60000}"
LOGDIR=paper/references/logs
mkdir -p "$LOGDIR" paper/references/results paper/quantum/results
python -c "import ase" 2>/dev/null || pip install -q ase

JOBLIST=(
  "paper/references/run_soc_benchmark.py I2-BDP"
  "paper/references/run_soc_benchmark.py I2-BDP-0"
  "paper/references/run_adiabatic.py BDP"
  "paper/references/run_adiabatic.py I2-BDP"
  "paper/quantum/thc_error.py BODIPY_core --ranks 22 33 44 --starts 8"
)

push() {
  git add paper/references/results/*.soc_benchmark.json paper/references/results/*.adiabatic* \
          paper/quantum/results/*.thc_error.json 2>/dev/null
  if git diff --cached --quiet; then return 0; fi
  git -c user.name="Parfait" -c user.email="26543855+pifparfait@users.noreply.github.com" \
    commit -q -m "Add SOC benchmark / adiabatic / THC-error results: $1 ($(date -u +%Y-%m-%d))"
  if ! git push -q origin "HEAD:refs/heads/results-tars-checks2"; then
    alt="results-tars-checks2-$(date -u +%Y%m%d-%H%M%S)"
    echo "Push rejected; not forcing. Publishing to $alt instead." >&2
    git push -q origin "HEAD:refs/heads/$alt"
  fi
}

run_one() {
  local tag; tag=$(echo "$1" | awk '{print $1"_"$2}' | sed 's|.*/||; s|\.py||')
  OMP_NUM_THREADS=$THREADS PSS_MAX_MEMORY_MB=$MEM_MB python $1 > "$LOGDIR/checks2_${tag}.log" 2>&1
  echo "finished $1 (exit $?)"
  ( flock 9; push "$tag" ) 9>"/tmp/psscreen_checks2_$USER.lock"   # parallel jobs push one at a time
}

for job in "${JOBLIST[@]}"; do
  while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do wait -n; done
  run_one "$job" &
  echo "started $job"
done
wait
push "all"
echo "all checks2 jobs finished"
