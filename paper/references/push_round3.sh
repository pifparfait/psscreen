#!/usr/bin/env bash
# Commit and push round-3 results to results-tars-round3 as Parfait. With --watch, repeat every hour while
# run_round3_tars.sh is running, then once more at the end.
# (run_round3_tars.sh passed all globs to one `git add`; a glob with no match made git add nothing.)
cd "$(dirname "$0")/../.."
push_once() {
  for f in paper/references/results/*.json paper/references/results/*.npz paper/references/results/*.xyz \
           paper/references/results/*.log validation/runs/set3/*.json validation/runs/set3/*.log \
           paper/references/logs/round3*.log; do
    [ -e "$f" ] && git add -f "$f"
  done
  git diff --cached --quiet && { echo "$(date) nothing new"; return 0; }
  git -c user.name="Parfait" -c user.email="26543855+pifparfait@users.noreply.github.com" \
    commit -q -m "Add round-3 review calculations ($(date -u +%Y-%m-%dT%H:%MZ))"
  git push -q origin "HEAD:refs/heads/results-tars-round3" && echo "$(date) pushed"
}
if [ "${1:-}" = "--watch" ]; then
  while pgrep -f run_round3_tars.sh >/dev/null; do push_once; sleep 3600; done
fi
push_once
