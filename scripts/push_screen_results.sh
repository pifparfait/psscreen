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
# Never force-push. If the remote branch has commits this clone does not have, publish to a new
# timestamped branch instead so that nothing already on GitHub is overwritten.
if ! git push -u origin "$BRANCH"; then
  alt="$BRANCH-$(date -u +%Y%m%d-%H%M%S)"
  echo "Push to $BRANCH rejected; not forcing. Publishing to $alt instead." >&2
  git push -u origin "HEAD:refs/heads/$alt"
  BRANCH="$alt"
fi
echo "Pushed screening results to branch $BRANCH."
