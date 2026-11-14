#!/usr/bin/env bash
# Screen a candidate CSV with one process per molecule (better throughput than one process with
# many threads: DFT/TDDFT in PySCF scales poorly beyond ~8 threads). Resumable.
# Usage: JOBS=8 THREADS=7 AUTO_PUSH=1 bash scripts/run_screen_parallel.sh screening/candidates_round1.csv [OUT_DIR]
set -euo pipefail
csv="${1:?usage: run_screen_parallel.sh CANDIDATES.csv [OUT_DIR]}"
cd "$(dirname "$0")/.."
out="${2:-screening/runs/$(basename "${csv%.csv}")}"
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate psscreen
export OMP_NUM_THREADS="${THREADS:-8}" OPENBLAS_NUM_THREADS="${THREADS:-8}" MKL_NUM_THREADS="${THREADS:-8}"
export CSV="$csv" OUT="$out" JOBS="${JOBS:-8}" AUTO_PUSH="${AUTO_PUSH:-0}" RESULTS_BRANCH="${RESULTS_BRANCH:-results-tars}"
mkdir -p "$out/inputs"
python - <<'PY'
import csv, os
rows = list(csv.DictReader(open(os.environ["CSV"])))
for r in rows:
    with open(os.path.join(os.environ["OUT"], "inputs", r["name"] + ".csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["name", "smiles"]); w.writeheader(); w.writerow(r)
PY
nohup bash -c '
  ls "$OUT"/inputs/*.csv | xargs -P "$JOBS" -I{} sh -c \
    "python screen.py {} --out \"$OUT\" --skip-existing > \"$OUT/\$(basename {} .csv).log\" 2>&1"
  python screen.py "$CSV" --out "$OUT" --skip-existing --summary "$OUT/ranking.md" > "$OUT/screen.log" 2>&1
  if [ "$AUTO_PUSH" = "1" ]; then bash scripts/push_screen_results.sh "$OUT"; fi
' > "$out/server_run.log" 2>&1 &
echo "started in background (PID $!); follow with: tail -f $out/*.log"
