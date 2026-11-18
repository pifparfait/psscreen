#!/usr/bin/env bash
# pi-controlled (12e,11o) active space (AVAS on the dipyrrin 2p_z) for the validation dyes on tars, then push the
# JSON results to the branch results-tars-pi as Parfait (never force-pushed). Runs JOBS dyes at a time and pushes
# after every finished dye, so partial results arrive early.
# Usage (inside a worktree of the paper branch, psscreen conda env active):
#   nohup bash paper/references/run_pi_tars.sh > paper/references/logs/pi.log 2>&1 &
set -uo pipefail
cd "$(dirname "$0")/../.."
THREADS="${THREADS:-8}"
JOBS="${JOBS:-4}"
MEM_MB="${PSS_MAX_MEMORY_MB:-40000}"
DYES=("$@")
[ ${#DYES[@]} -eq 0 ] && DYES=(BDP I2-BDP BODIPY_core I-BDP I3-BDP azaBDP_Ph4 BDP-0 I-BDP-0 I2-BDP-0)
LOGDIR=paper/references/logs
mkdir -p "$LOGDIR" paper/references/results

push() {
  git add paper/references/results/*.pi_cas_12_11.json
  if git diff --cached --quiet; then return 0; fi
  git -c user.name="Parfait" -c user.email="26543855+pifparfait@users.noreply.github.com" \
    commit -q -m "Add pi-controlled CAS(12,11) results: $1 ($(date -u +%Y-%m-%d))"
  if ! git push -q origin "HEAD:refs/heads/results-tars-pi"; then
    alt="results-tars-pi-$(date -u +%Y%m%d-%H%M%S)"
    echo "Push rejected; not forcing. Publishing to $alt instead." >&2
    git push -q origin "HEAD:refs/heads/$alt"
  fi
}

run_one() {
  OMP_NUM_THREADS=$THREADS PSS_MAX_MEMORY_MB=$MEM_MB \
    python paper/references/run_pi_cas.py "$1" > "$LOGDIR/${1}_pi.log" 2>&1
  echo "finished $1 (exit $?)"
}

for d in "${DYES[@]}"; do
  while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do wait -n; push "progress"; done
  run_one "$d" &
  echo "started $d"
done
wait
push "${DYES[*]}"
echo "all pi runs finished"
