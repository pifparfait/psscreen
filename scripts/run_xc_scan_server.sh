#!/usr/bin/env bash
# Functional scan (B3LYP baseline already exists): PBE0 and BHandHLYP on the 9 validation dyes.
# One process per dye and functional, THREADS threads each; resumable. Usage:
#   JOBS=6 THREADS=8 AUTO_PUSH=1 bash scripts/run_xc_scan_server.sh
set -euo pipefail
cd "$(dirname "$0")/../validation"
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate psscreen
export OMP_NUM_THREADS="${THREADS:-8}" PSS_EWIN_SINGLETS="${PSS_EWIN_SINGLETS:-60}"
export PSS_MAX_MEMORY_MB="${PSS_MAX_MEMORY_MB:-30000}"
JOBS="${JOBS:-6}"; AUTO_PUSH="${AUTO_PUSH:-0}"
mkdir -p runs_xc
nohup bash -c '
  for xc in bhandhlyp pbe0; do for n in I2-BDP-0 BDP-0 I-BDP-0 I2-BDP I3-BDP I-BDP BDP BODIPY_core azaBDP_Ph4; do echo "$xc $n"; done; done |
    xargs -P '"$JOBS"' -L1 bash -c '"'"'python run_xc_scan.py $0 $1 >> runs_xc/scan_$0_$1.log 2>&1'"'"'
  python report_xc.py > ../docs/XC_SCAN.md 2>&1
  if [ "'"$AUTO_PUSH"'" = "1" ]; then bash ../scripts/push_results.sh; fi
' > runs_xc/server_run.log 2>&1 &
echo "started (PID $!); follow with: tail -f validation/runs_xc/scan_*.log"
