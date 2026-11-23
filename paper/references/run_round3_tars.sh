#!/usr/bin/env bash
# Third mock-review round, calculations on tars:
#   S0 Hessians of BDP and I2-BDP (multimode Franck-Condon density of states),
#   effective SOC operator on the all-electron TDA states of the SOMF benchmark (I2-BDP, I2-BDP-0),
#   TDA gradients of S1 and T1-T3 at the S0 minimum (vertical-gradient relaxation, multimode line shapes),
#   S0/T1/S1 optimizations, S0 Hessian and TDA gradients of aza-BDP-Br2 (relaxation of the near-resonant S1-T2
#   channel; T2 itself cannot be held as a DeltaSCF determinant),
#   AVAS pi-space SA-CASSCF/SC-NEVPT2 of BDP and I2-BDP at the DFT S0 minimum instead of the xTB geometry,
#   the missing c1 conformers of BDP-0 and I-BDP-0.
# Results are pushed to the branch results-tars-round3 as Parfait after every finished job (never force-pushed).
# Usage (inside ~/psscreen-paper with the paper branch checked out, psscreen env active):
#   mkdir -p paper/references/logs && nohup bash paper/references/run_round3_tars.sh > paper/references/logs/round3.log 2>&1 &
set -uo pipefail
export PYTHONUNBUFFERED=1
cd "$(dirname "$0")/../.."
THREADS="${THREADS:-16}"
JOBS="${JOBS:-4}"
export PSS_MAX_MEMORY_MB="${PSS_MAX_MEMORY_MB:-50000}"
LOGDIR=paper/references/logs
mkdir -p "$LOGDIR" paper/references/results
python -c "import ase" 2>/dev/null || pip install -q ase

JOBLIST=(
  "chain_azaBr2|python paper/references/run_adiabatic.py azaBDP_Br2 --states S0 T1 S1 && python paper/references/run_s0_hessian.py azaBDP_Br2 && python paper/references/run_vertical_gradients.py azaBDP_Br2 --states S1 T1 T2"
  "hess_BDP|python paper/references/run_s0_hessian.py BDP"
  "hess_I2-BDP|python paper/references/run_s0_hessian.py I2-BDP"
  "vg_BDP|python paper/references/run_vertical_gradients.py BDP --states S1 T1"
  "vg_I2-BDP|python paper/references/run_vertical_gradients.py I2-BDP --states S1 T1 T2"
  "soc_I2-BDP|python paper/references/run_soc_benchmark.py I2-BDP --effective"
  "soc_I2-BDP-0|python paper/references/run_soc_benchmark.py I2-BDP-0 --effective"
  "pi_BDP_dft|python paper/references/run_pi_cas.py BDP --geom paper/references/results/BDP.adiabatic_S0.xyz --tag _dftS0 --skip-avas-block"
  "pi_I2-BDP_dft|python paper/references/run_pi_cas.py I2-BDP --geom paper/references/results/I2-BDP.adiabatic_S0.xyz --tag _dftS0 --skip-avas-block"
  "conf_c1|python screen.py validation/set3_c1.csv --out validation/runs/set3 --skip-existing"
)

push() {
  git add -f paper/references/results/*.json paper/references/results/*.npz paper/references/results/*.xyz \
            paper/references/results/*.log validation/runs/set3/*.json validation/runs/set3/*.log 2>/dev/null
  if git diff --cached --quiet; then return 0; fi
  git -c user.name="Parfait" -c user.email="26543855+pifparfait@users.noreply.github.com" \
    commit -q -m "Add round-3 review calculations: $1 ($(date -u +%Y-%m-%d))"
  if ! git push -q origin "HEAD:refs/heads/results-tars-round3"; then
    alt="results-tars-round3-$(date -u +%Y%m%d-%H%M%S)"
    echo "Push rejected; not forcing. Publishing to $alt instead." >&2
    git push -q origin "HEAD:refs/heads/$alt"
  fi
}

run_one() {
  local tag="$1" cmd="$2"
  OMP_NUM_THREADS=$THREADS bash -c "$cmd" > "$LOGDIR/round3_${tag}.log" 2>&1
  echo "finished $tag (exit $?)"
  ( flock 9; push "$tag" ) 9>"/tmp/psscreen_round3_$USER.lock"
}

for entry in "${JOBLIST[@]}"; do
  tag="${entry%%|*}"; cmd="${entry#*|}"
  if [ "$#" -gt 0 ] && [[ " $* " != *" $tag "* ]]; then continue; fi
  while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do wait -n; done
  run_one "$tag" "$cmd" &
  echo "started $tag"
done
wait
push "all"
echo "all round-3 jobs finished"
