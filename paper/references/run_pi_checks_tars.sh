#!/usr/bin/env bash
# Checks of the pi-controlled CAS(12,11) results for BDP and I2-BDP on tars: def2-TZVP basis, and state averaging
# over singlets + triplets. Pushes JSON results to results-tars-pi-checks as Parfait (never force-pushed), after
# every finished job.
# Usage (inside a worktree of the paper branch, psscreen conda env active):
#   nohup bash paper/references/run_pi_checks_tars.sh > paper/references/logs/pi_checks.log 2>&1 &
set -uo pipefail
cd "$(dirname "$0")/../.."
THREADS="${THREADS:-12}"
JOBS="${JOBS:-2}"
MEM_MB="${PSS_MAX_MEMORY_MB:-60000}"
LOGDIR=paper/references/logs
mkdir -p "$LOGDIR" paper/references/results
JOBLIST=("BDP --basis def2-tzvp" "I2-BDP --basis def2-tzvp" "BDP --sa-mix" "I2-BDP --sa-mix")

push() {
  git add paper/references/results/*.pi_cas_12_11_*.json 2>/dev/null
  if git diff --cached --quiet; then return 0; fi
  git -c user.name="Parfait" -c user.email="26543855+pifparfait@users.noreply.github.com" \
    commit -q -m "Add pi-CAS checks (TZVP / SA singlets+triplets): $1 ($(date -u +%Y-%m-%d))"
  if ! git push -q origin "HEAD:refs/heads/results-tars-pi-checks"; then
    alt="results-tars-pi-checks-$(date -u +%Y%m%d-%H%M%S)"
    echo "Push rejected; not forcing. Publishing to $alt instead." >&2
    git push -q origin "HEAD:refs/heads/$alt"
  fi
}

run_one() {
  local tag; tag=$(echo "$1" | tr ' ' '_' | tr -d '-')
  OMP_NUM_THREADS=$THREADS PSS_MAX_MEMORY_MB=$MEM_MB \
    python paper/references/run_pi_cas.py $1 --skip-avas-block > "$LOGDIR/pi_check_${tag}.log" 2>&1
  echo "finished $1 (exit $?)"
}

for j in "${JOBLIST[@]}"; do
  while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do wait -n; push "progress"; done
  run_one "$j" &
  echo "started $j"
done
wait
push "all"
echo "all pi checks finished"
