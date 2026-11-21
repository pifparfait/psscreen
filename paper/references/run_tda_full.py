"""Convergence of the pipeline triplet energies with the TDA orbital window.

The pipeline computes triplets with an explicit TDA in an orbital window (10 eV around the frontier orbitals,
at most 40 occupied and 60 virtual orbitals, XC kernel on PySCF grid level 1). The upper triplets T2, T3 decide
the ISC rate of azaBDP_Br2 and the discussion of I2-BDP-0, so here they are recomputed for the same geometry,
functional and basis (B3LYP/def2-SVP, density fitting, def2 ECP on I/Br... as in the pipeline) with
  (1) the pipeline window, recomputed in this run (control),
  (2) a 20 eV window without orbital caps,
  (3) the full TDA of PySCF (Davidson, all orbital pairs, kernel on grid level 2).
Usage (repository root): python paper/references/run_tda_full.py NAME [--nstates 6]
Writes paper/references/results/NAME.tda_full.json.
"""
import argparse
import json
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from pyscf import dft, gto, tdscf  # noqa: E402

from psscreen.backends import _ecp  # noqa: E402
from psscreen.geometry import read_xyz, to_pyscf_atom  # noqa: E402
from psscreen.tdwin import window_td  # noqa: E402

EV = 27.211386245988
ap = argparse.ArgumentParser()
ap.add_argument("name")
ap.add_argument("--nstates", type=int, default=6)
a = ap.parse_args()

CANDIDATES = [os.path.join(ROOT, "validation", "runs", f"{a.name}"),
              os.path.join(ROOT, "validation", "runs", "set3", f"{a.name}_c0"),
              os.path.join(ROOT, "screening", "runs", "candidates_round1", f"{a.name}")]
base = next(p for p in CANDIDATES if os.path.exists(p + ".xyz"))
pipeline = json.load(open(base + ".json")) if os.path.exists(base + ".json") else {}

sym, xyz = read_xyz(base + ".xyz")
atom = to_pyscf_atom(sym, xyz)
mem = int(os.environ.get("PSS_MAX_MEMORY_MB", "60000"))
mol = gto.M(atom=atom, basis="def2-svp", ecp=_ecp(atom, "def2-svp"), verbose=4, max_memory=mem)
mf = dft.RKS(mol, xc="b3lyp").density_fit()
mf.grids.level = 2
mf.conv_tol = 1e-9
mf.kernel()
assert mf.converged
out = {"name": a.name, "geometry": base + ".xyz", "nao": mol.nao,
       "nocc": mol.nelectron // 2, "pipeline_E_T_eV": pipeline.get("E_triplets_eV"),
       "pipeline_window": {k: pipeline.get("meta", {}).get(k) for k in ("ewin_eV", "n_occ_win", "n_vir_win")}}

t0 = time.time()
r = window_td(mf, ewin_ev=10.0, nstates=a.nstates, rpa_singlets=False, max_occ=40, max_vir=60, grid_level=1)
out["window10"] = {"E_T_eV": [float(e) * EV for e in r["e_T"]], "time_s": time.time() - t0}
print("window 10 eV:", out["window10"]["E_T_eV"], flush=True)

t0 = time.time()
r = window_td(mf, ewin_ev=20.0, nstates=a.nstates, rpa_singlets=False, max_occ=999, max_vir=999, grid_level=1)
out["window20"] = {"E_T_eV": [float(e) * EV for e in r["e_T"]], "time_s": time.time() - t0}
print("window 20 eV:", out["window20"]["E_T_eV"], flush=True)

t0 = time.time()
td = tdscf.TDA(mf)
td.singlet = False
td.nstates = a.nstates + 2
td.conv_tol = 1e-6
td.kernel()
out["full"] = {"E_T_eV": [float(e) * EV for e in td.e[: a.nstates]], "converged": [bool(c) for c in td.converged],
               "time_s": time.time() - t0}
print("full TDA:", out["full"]["E_T_eV"], flush=True)

os.makedirs(os.path.join(ROOT, "paper", "references", "results"), exist_ok=True)
with open(os.path.join(ROOT, "paper", "references", "results", f"{a.name}.tda_full.json"), "w") as f:
    json.dump(out, f, indent=1)
