#!/usr/bin/env bash
# Commit and push server validation results to the `results-tars` branch.
# Requires: git remote "origin" with write access (e.g. a deploy key) and
# git user.name / user.email configured.
set -euo pipefail
cd "$(dirname "$0")/.."

BRANCH="${RESULTS_BRANCH:-results-tars}"

if [ ! -d .git ]; then
  echo "Not a git repository. Run: git init -b main && git remote add origin <url>" >&2
  exit 1
fi

git checkout -B "$BRANCH"
git add validation/runs/*.json docs/VALIDATION_server.md 2>/dev/null || true
if git diff --cached --quiet; then
  echo "Nothing new to commit."
else
  git commit -m "Add server validation results ($(date -u +%Y-%m-%d))"
fi
git push -u origin "$BRANCH" --force-with-lease || git push -u origin "$BRANCH" --force
echo "Pushed results to branch $BRANCH."
