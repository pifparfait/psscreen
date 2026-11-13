#!/usr/bin/env bash
# Run the whole validation (classical + quantum backend + CASCI/NEVPT2 reference) in the
# background on a server, logging to validation/runs/*.log. Resumable: finished molecules are
# skipped. Usage:
#   bash scripts/run_validation_server.sh            # uses all cores, full singlet space
#   OMP_NUM_THREADS=8 PSS_EWIN_SINGLETS=25 bash scripts/run_validation_server.sh   # low memory
# Environment:
#   FORCE=1       recompute everything, ignoring result JSONs already present
#   AUTO_PUSH=1   when finished, run scripts/push_results.sh (needs a git remote with write access)
set -euo pipefail
cd "$(dirname "$0")/../validation"
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate psscreen
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-$(nproc)}"
export PSS_EWIN_SINGLETS="${PSS_EWIN_SINGLETS:-60}"
export PSS_MAX_MEMORY_MB="${PSS_MAX_MEMORY_MB:-100000}"  # PySCF memory limit (MB)
export FORCE="${FORCE:-0}" AUTO_PUSH="${AUTO_PUSH:-0}"
mkdir -p runs
if [ "$FORCE" = "1" ]; then
  rm -f runs/*.json
  git checkout -- runs/*.xyz 2>/dev/null || true
fi
nohup bash -c '
  for m in BODIPY_core BDP I2-BDP azaBDP_Ph4; do
    [ -f runs/$m.json ] || python run_validation.py $m >> runs/validation.log 2>&1
  done
  for m in BODIPY_core BDP I2-BDP azaBDP_Ph4; do
    [ -f runs/$m.quantum.json ] || python run_quantum.py $m 4 6 > runs/quantum_$m.log 2>&1
  done
  for m in BODIPY_core BDP I2-BDP azaBDP_Ph4; do
    [ -f runs/$m.casci_ref.json ] || python run_casci_reference.py $m --nevpt2 > runs/casci_$m.log 2>&1
  done
  python report.py > ../docs/VALIDATION_server.md 2>&1
  if [ "$AUTO_PUSH" = "1" ]; then bash ../scripts/push_results.sh; fi
' > runs/server_run.log 2>&1 &
echo "started in background (PID $!); follow with: tail -f validation/runs/validation.log"
