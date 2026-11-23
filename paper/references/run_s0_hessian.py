"""Ground-state Hessian at the DFT S0 minimum, for a multimode Franck-Condon weighted density of states.

The MLJ rate of the paper uses an assumed single-mode line shape (lambda_s = 0.10 eV, hbar*omega = 0.17 eV, S = 0.8).
This script computes the B3LYP/def2-SVP (density fitting, def2 ECP on iodine, grid level 3) Hessian at the S0
minimum of run_adiabatic.py (paper/references/results/NAME.adiabatic_S0.xyz). Analytic second derivatives are used
when PySCF provides them for this setup; otherwise the Hessian is built by central finite differences of analytic
gradients (step 0.005 Bohr), check-pointed row by row. The S1 and T1 minima of run_adiabatic.py are then projected on
the normal modes by paper/analysis/fcwd_multimode.py (displaced-oscillator model, S0 frequencies).

Usage (repository root): PYTHONUNBUFFERED=1 python paper/references/run_s0_hessian.py NAME [--fd]
Environment: OMP_NUM_THREADS, PSS_MAX_MEMORY_MB (default 40000).
Writes paper/references/results/NAME.s0_hessian.npz (Hessian in Ha/Bohr^2, geometry in Angstrom) and .json.
"""
import argparse
import json
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
import numpy as np  # noqa: E402
from pyscf import dft, gto  # noqa: E402

from psscreen.backends import _ecp  # noqa: E402
from psscreen.geometry import read_xyz, to_pyscf_atom  # noqa: E402

BOHR_A = 0.52917721092
ap = argparse.ArgumentParser()
ap.add_argument("name")
ap.add_argument("--fd", action="store_true", help="force finite differences of gradients")
ap.add_argument("--step", type=float, default=0.005, help="Bohr")
a = ap.parse_args()
RES = os.path.join(ROOT, "paper", "references", "results")
NPZ = os.path.join(RES, f"{a.name}.s0_hessian.npz")
CKPT = os.path.join(RES, f"{a.name}.s0_hessian_fd_partial.npz")
MEM = int(os.environ.get("PSS_MAX_MEMORY_MB", 40000))
sym, xyz = read_xyz(os.path.join(RES, f"{a.name}.adiabatic_S0.xyz"))
sym = [s.capitalize() for s in sym]
xyz = np.asarray(xyz, float)
n = len(sym)


def mf_at(x):
    atom = to_pyscf_atom(sym, x)
    mol = gto.M(atom=atom, basis="def2-svp", ecp=_ecp(atom, "def2-svp"), verbose=0, max_memory=MEM)
    mf = dft.RKS(mol, xc="b3lyp").density_fit()
    mf.grids.level = 3
    mf.conv_tol = 1e-11
    mf.conv_tol_grad = 1e-7
    return mf


t0 = time.time()
mf = mf_at(xyz)
e0 = mf.kernel()
g0 = mf.nuc_grad_method().kernel()
method = None
H = None
if not a.fd:
    try:
        h4 = mf.Hessian().kernel()                       # (n, n, 3, 3)
        H = h4.transpose(0, 2, 1, 3).reshape(3 * n, 3 * n)
        method = "analytic"
    except Exception as exc:                             # noqa: BLE001
        print("analytic Hessian unavailable:", repr(exc), flush=True)
if H is None:
    method = f"finite differences of analytic gradients, step {a.step} Bohr"
    H = np.full((3 * n, 3 * n), np.nan)
    if os.path.exists(CKPT):
        H = np.load(CKPT)["H"]
    dm0 = mf.make_rdm1()
    for k in range(3 * n):
        if not np.isnan(H[k]).any():
            continue
        rows = []
        for sgn in (1, -1):
            x = xyz.copy()
            x[k // 3, k % 3] += sgn * a.step * BOHR_A
            m = mf_at(x)
            m.kernel(dm0=dm0)
            rows.append(np.asarray(m.nuc_grad_method().kernel()).ravel())
        H[k] = (rows[0] - rows[1]) / (2 * a.step)
        np.savez(CKPT, H=H)
        print(f"row {k + 1}/{3 * n} done, {time.time() - t0:.0f} s", flush=True)
    H = 0.5 * (H + H.T)
np.savez(NPZ, H=H, xyz_A=xyz, symbols=np.array(sym), grad=np.asarray(g0))
json.dump({"name": a.name, "method": "B3LYP/def2-SVP DF, grid 3, " + method, "E_S0_Ha": float(e0),
           "max_abs_grad_Ha_bohr": float(np.abs(g0).max()), "time_s": time.time() - t0},
          open(os.path.join(RES, f"{a.name}.s0_hessian.json"), "w"), indent=1)
print(a.name, "done", method, round(time.time() - t0), "s")
