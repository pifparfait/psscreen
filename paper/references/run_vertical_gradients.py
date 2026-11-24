"""TDA excited-state gradients at the DFT S0 minimum (vertical-gradient model of the excited-state relaxation).

For the near-resonant S1-T2 channel of aza-BDP-Br2 the second triplet cannot be optimized as a DeltaSCF determinant
(the MOM solution collapses towards T1). In the vertical-gradient (displaced-oscillator) model the relaxation of a
state X from the S0 minimum is set by its gradient there and the S0 Hessian H:
    d_X = -H^-1 g_X  (translations and rotations projected out),   lambda_X = (1/2) g_X^T H^-1 g_X,
and the displacement between two states, d_S1 - d_Tn projected on the mass-weighted normal modes, gives the
Huang-Rhys factors of that S1 -> T_n channel (paper/analysis/fcwd_multimode.py).
This script computes, at paper/references/results/NAME.adiabatic_S0.xyz, the B3LYP/def2-SVP (no density fitting:
PySCF has no density-fitted TDDFT gradients; def2 ECP on iodine, grid level 2 as in the pipeline) TDA gradients
of S1 and of T1, T2, T3 (singlet=False), or of the states requested, with the TDA excitation energies.

Usage (repository root): PYTHONUNBUFFERED=1 python paper/references/run_vertical_gradients.py NAME [--geom XYZ]
                         [--states S1 T1 T2 T3]
Environment: OMP_NUM_THREADS, PSS_MAX_MEMORY_MB (default 40000).
Writes paper/references/results/NAME.vertical_gradients.npz (gradients in Ha/Bohr) and .json (energies, timings).
"""
import argparse
import json
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
import numpy as np  # noqa: E402
from pyscf import dft, gto, tdscf  # noqa: E402

from psscreen.backends import _ecp  # noqa: E402
from psscreen.geometry import read_xyz, to_pyscf_atom  # noqa: E402

EV = 27.211386245988
ap = argparse.ArgumentParser()
ap.add_argument("name")
ap.add_argument("--geom", default=None)
ap.add_argument("--states", nargs="+", default=["S1", "T1", "T2", "T3"], help="S1 and/or T1, T2, T3, ...")
a = ap.parse_args()
RES = os.path.join(ROOT, "paper", "references", "results")
MEM = int(os.environ.get("PSS_MAX_MEMORY_MB", 40000))
geo = a.geom or os.path.join(RES, f"{a.name}.adiabatic_S0.xyz")
sym, xyz = read_xyz(geo)
sym = [s.capitalize() for s in sym]
atom = to_pyscf_atom(sym, np.asarray(xyz, float))
mol = gto.M(atom=atom, basis="def2-svp", ecp=_ecp(atom, "def2-svp"), verbose=0, max_memory=MEM)
t0 = time.time()
mf = dft.RKS(mol, xc="b3lyp")
mf.grids.level = 2
mf.conv_tol = 1e-10
mf.kernel()
assert mf.converged
out = {"name": a.name, "geometry": os.path.relpath(geo, ROOT), "method": "B3LYP/def2-SVP, no DF, grid 2, TDA",
       "E_S0_Ha": float(mf.e_tot), "timings_s": {"scf": time.time() - t0}}
grads = {"S0": mf.nuc_grad_method().kernel(), "xyz_A": np.asarray(xyz, float)}
out["timings_s"]["S0_grad"] = time.time() - t0

want = {"S": sorted(int(x[1:]) - 1 for x in a.states if x[0] == "S"),
        "T": sorted(int(x[1:]) - 1 for x in a.states if x[0] == "T")}
for label, singlet in (("S", True), ("T", False)):
    roots = want[label]
    if not roots:
        continue
    # Only the roots whose gradients we need must converge; the extra guard root
    # often stalls at tight tolerance, so retry with a larger subspace if needed.
    nst = max(roots) + 3
    for tol, cyc, space in ((1e-7, 200, 40), (1e-6, 400, 80)):
        td = tdscf.TDA(mf)
        td.singlet, td.nstates, td.conv_tol, td.max_cycle, td.max_space = singlet, nst, tol, cyc, space
        td.kernel()
        if all(td.converged[r] for r in roots):
            break
    assert all(td.converged[r] for r in roots), f"TDA ({label}) roots {roots} not converged"
    out[f"tda_{label}_conv_tol"] = tol
    out[f"E_{label}_tda_eV"] = (np.asarray(td.e) * EV).tolist()
    for r in roots:
        t1 = time.time()
        g = td.nuc_grad_method().kernel(state=r + 1)
        grads[f"{label}{r + 1}"] = g
        out["timings_s"][f"{label}{r + 1}_grad"] = time.time() - t1
        print(label, r + 1, f"E = {td.e[r] * EV:.4f} eV, |g| = {np.linalg.norm(g):.4e} Ha/Bohr", flush=True)
        np.savez(os.path.join(RES, f"{a.name}.vertical_gradients.npz"), symbols=np.array(sym), **grads)
        json.dump(out, open(os.path.join(RES, f"{a.name}.vertical_gradients.json"), "w"), indent=1)
out["timings_s"]["total"] = time.time() - t0
json.dump(out, open(os.path.join(RES, f"{a.name}.vertical_gradients.json"), "w"), indent=1)
print(a.name, "done", round(time.time() - t0), "s")
