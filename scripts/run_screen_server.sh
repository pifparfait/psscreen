#!/usr/bin/env bash
# Screen a candidate CSV in the background on a server (resumable), then optionally push results.
# Usage: bash scripts/run_screen_server.sh screening/candidates_round1.csv [OUT_DIR]
# Environment: OMP_NUM_THREADS (default: all cores), AUTO_PUSH=1 (push to RESULTS_BRANCH when done)
set -euo pipefail
csv="${1:?usage: run_screen_server.sh CANDIDATES.csv [OUT_DIR]}"
cd "$(dirname "$0")/.."
out="${2:-screening/runs/$(basename "${csv%.csv}")}"
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate psscreen
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-$(nproc)}"
export CSV="$csv" OUT="$out" AUTO_PUSH="${AUTO_PUSH:-0}"
export RESULTS_BRANCH="${RESULTS_BRANCH:-results-tars}"
mkdir -p "$out"
nohup bash -c '
  python screen.py "$CSV" --out "$OUT" --skip-existing --summary "$OUT/ranking.md" > "$OUT/screen.log" 2>&1
  if [ "$AUTO_PUSH" = "1" ]; then bash scripts/push_screen_results.sh "$OUT"; fi
' > "$out/server_run.log" 2>&1 &
echo "started in background (PID $!); follow with: tail -f $out/screen.log"
