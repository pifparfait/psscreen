#!/usr/bin/env bash
# Commit screening results (JSON, ranking, geometries) from OUT_DIR and push to RESULTS_BRANCH.
set -euo pipefail
cd "$(dirname "$0")/.."
out="${1:?usage: push_screen_results.sh OUT_DIR}"
BRANCH="${RESULTS_BRANCH:-results-tars}"
[ -d .git ] || { echo "Not a git repository" >&2; exit 1; }
git checkout -B "$BRANCH"
git add "$out"/*.json "$out"/*.xyz "$out"/ranking.md "$out"/*.log 2>/dev/null || true
if git diff --cached --quiet; then echo "Nothing new to commit."
else git commit -m "Add screening results for $(basename "$out") ($(date -u +%Y-%m-%d))"; fi
git push -u origin "$BRANCH" --force-with-lease || git push -u origin "$BRANCH" --force
echo "Pushed screening results to branch $BRANCH."
