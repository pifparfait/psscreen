#!/usr/bin/env bash
# Run the SA-CASSCF/SC-NEVPT2 singlet+triplet references for all nine validation dyes in parallel,
# then push paper/references/results/*.json to the branch `results-tars-paper` (never force-pushed).
#
# Usage (inside a worktree of the paper branch, with the psscreen conda env active):
#   bash paper/references/run_all_tars.sh            # CAS(12,11), 9 jobs x 7 threads
#   CAS=14,13 bash paper/references/run_all_tars.sh BDP-0 I-BDP-0 I2-BDP-0
set -uo pipefail
cd "$(dirname "$0")/../.."

CAS="${CAS:-12,11}"
THREADS="${THREADS:-7}"
MEM_MB="${PSS_MAX_MEMORY_MB:-24000}"
DYES=("$@")
[ ${#DYES[@]} -eq 0 ] && DYES=(BODIPY_core BDP I-BDP I2-BDP I3-BDP azaBDP_Ph4 BDP-0 I-BDP-0 I2-BDP-0)
LOGDIR=paper/references/logs
mkdir -p "$LOGDIR" paper/references/results

pids=()
for d in "${DYES[@]}"; do
  OMP_NUM_THREADS=$THREADS PSS_MAX_MEMORY_MB=$MEM_MB \
    python paper/references/run_sa_nevpt2.py "$d" --cas "$CAS" --nsing 3 --ntrip 3 \
    > "$LOGDIR/${d}_${CAS/,/_}.log" 2>&1 &
  pids+=($!)
  echo "started $d (pid $!)"
done
for p in "${pids[@]}"; do wait "$p"; done
echo "all jobs finished"

# Push the results as Parfait to a dedicated branch.
git add paper/references/results/*.json
if git diff --cached --quiet; then echo "Nothing new to commit."; exit 0; fi
git -c user.name="Parfait" -c user.email="26543855+pifparfait@users.noreply.github.com" \
  commit -m "Add SA-CASSCF/SC-NEVPT2 singlet and triplet references, CAS($CAS) ($(date -u +%Y-%m-%d))"
if ! git push origin "HEAD:refs/heads/results-tars-paper"; then
  alt="results-tars-paper-$(date -u +%Y%m%d-%H%M%S)"
  echo "Push rejected; not forcing. Publishing to $alt instead." >&2
  git push origin "HEAD:refs/heads/$alt"
fi
