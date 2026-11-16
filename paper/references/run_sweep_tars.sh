#!/usr/bin/env bash
# Active-space convergence sweep for BDP and I2-BDP on tars, then push the JSON results to the branch
# results-tars-sweep as Parfait (never force-pushed). Runs the two dyes in parallel.
# Usage (inside a worktree of the paper branch, psscreen conda env active):
#   nohup bash paper/references/run_sweep_tars.sh > paper/references/logs/sweep.log 2>&1 &
set -uo pipefail
cd "$(dirname "$0")/../.."
THREADS="${THREADS:-14}"
MEM_MB="${PSS_MAX_MEMORY_MB:-60000}"
SIZES="${SIZES:-8,10,12,14,16}"
DYES=("$@")
[ ${#DYES[@]} -eq 0 ] && DYES=(BDP I2-BDP)
LOGDIR=paper/references/logs
mkdir -p "$LOGDIR" paper/references/results
pids=()
for d in "${DYES[@]}"; do
  OMP_NUM_THREADS=$THREADS PSS_MAX_MEMORY_MB=$MEM_MB \
    python paper/references/run_cas_sweep.py "$d" --sizes "$SIZES" > "$LOGDIR/${d}_sweep.log" 2>&1 &
  pids+=($!)
  echo "started $d (pid $!)"
done
for p in "${pids[@]}"; do wait "$p"; done
echo "all sweeps finished"
git add paper/references/results/*.cas_sweep.json
if git diff --cached --quiet; then echo "Nothing new to commit."; exit 0; fi
git -c user.name="Parfait" -c user.email="26543855+pifparfait@users.noreply.github.com" \
  commit -m "Add active-space convergence sweep for ${DYES[*]} ($(date -u +%Y-%m-%d))"
if ! git push origin "HEAD:refs/heads/results-tars-sweep"; then
  alt="results-tars-sweep-$(date -u +%Y%m%d-%H%M%S)"
  echo "Push rejected; not forcing. Publishing to $alt instead." >&2
  git push origin "HEAD:refs/heads/$alt"
fi
