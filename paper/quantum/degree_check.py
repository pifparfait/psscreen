"""Check of the filter-degree model used in decision_cost.py against our own Chebyshev filters.

decision_cost.py uses d = c (lambda'/Delta) sqrt(ln(1/eps)) + 321, with c taken from the preprint's empirical
step-function fit (d = 4.7571 lambda'/Delta + 321 at eps = 0.01). Here we build the actual filter
sqrt(rho_MLJ/rho_max) of psscreen.isc_filtered.filter_poly for a range of block-encoding normalizations alpha and
find the smallest even degree that reaches a uniform error eps, then fit d = a (alpha/Delta) + b.

It also computes, for the BODIPY-core CAS(12,11) Hamiltonian, the identity coefficient c0 = Tr H / 2^n of the
qubit Hamiltonian and the CASCI ground-state energy E0, so that the normalization of the shifted operator
(H - E_c)/lambda' with E_c just below E0 is lambda' = lambda + |E0 - c0| (lambda excludes the identity term).

Usage (from the repository root): python paper/quantum/degree_check.py
Writes paper/quantum/results/degree_check.json.
"""
import json
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
import numpy as np  # noqa: E402
from numpy.polynomial import chebyshev as C  # noqa: E402

from psscreen.isc_filtered import _rho_vec  # noqa: E402
from psscreen.quantum_pl import HARTREE_EV  # noqa: E402
from psscreen.rates import KB_EV, MLJParams  # noqa: E402

P = MLJParams()
DELTA_HA = math.sqrt(2 * P.lambda_s * KB_EV * P.T) / HARTREE_EV
GRID = np.linspace(-1.5, 0.5, 2001)
RHO_MAX = float(np.max(_rho_vec(GRID, P)))
GAP_HA = 2.5 / HARTREE_EV            # E_S1 - E_c, a typical S1 energy above the ground state


def f_factory(alpha):
    def f(x):
        dE = (alpha * np.abs(np.asarray(x, dtype=float)) - GAP_HA) * HARTREE_EV
        return np.sqrt(np.clip(_rho_vec(dE.ravel(), P).reshape(np.shape(dE)), 0, None) / RHO_MAX)
    return f


def err_at(f, deg, xs, fx):
    c = C.chebinterpolate(f, deg)
    return float(np.max(np.abs(C.chebval(xs, c) - fx)))


def min_degree(alpha, eps):
    f = f_factory(alpha)
    x0 = GAP_HA / alpha
    w = 40 * DELTA_HA / alpha
    xs = np.concatenate([np.linspace(-1, 1, 20001), np.linspace(x0 - w, x0 + 3 * w, 40001)])
    xs = np.clip(xs, -1, 1)
    fx = f(xs)
    lo, hi = 2, 64
    while err_at(f, hi, xs, fx) > eps:
        lo, hi = hi, 2 * hi
    while hi - lo > 4:
        mid = (lo + hi) // 2
        mid += mid % 2
        if err_at(f, mid, xs, fx) <= eps:
            hi = mid
        else:
            lo = mid
    return hi


def identity_and_e0():
    from pyscf import dft, gto, mcscf
    from psscreen.backends import _ecp
    from psscreen.geometry import read_xyz, to_pyscf_atom
    sym, xyz = read_xyz(os.path.join(ROOT, "validation", "runs", "BODIPY_core.xyz"))
    atom = to_pyscf_atom(sym, xyz)
    mol = gto.M(atom=atom, basis="def2-svp", ecp=_ecp(atom, "def2-svp"), verbose=0)
    mf = dft.RKS(mol, xc="b3lyp").density_fit()
    mf.grids.level = 2
    mf.kernel()
    mc = mcscf.CASCI(mf, 11, 12)
    h1, ecore = mc.get_h1eff()
    from pyscf import ao2mo
    g = ao2mo.restore(1, mc.get_h2eff(), 11)
    # Tr H / 2^n for H = sum h_pq E_pq + 1/2 sum (pq|rs) a+_p a+_r a_s a_q over spin-orbitals
    c0 = ecore + np.trace(h1) + 0.5 * np.einsum("ppqq->", g) - 0.25 * np.einsum("pqqp->", g)
    mc.fcisolver.nroots = 1
    e0 = mc.kernel()[0]
    return float(c0), float(e0), float(ecore)


def main():
    out = {"Delta_eV": DELTA_HA * HARTREE_EV, "gap_eV": GAP_HA * HARTREE_EV, "fits": {}}
    ratios = [100, 200, 400, 800, 1600]
    for eps in (1e-2, 2e-3, 1e-4):
        ds = []
        for r in ratios:
            d = min_degree(r * DELTA_HA, eps)
            ds.append(d)
            print(f"eps={eps:.0e} alpha/Delta={r:5d} degree={d}", flush=True)
        a, b = np.polyfit(ratios, ds, 1)
        out["fits"][f"{eps:.0e}"] = {"alpha_over_Delta": ratios, "degree": ds, "slope": float(a),
                                      "intercept": float(b),
                                      "slope_over_sqrt_ln": float(a / math.sqrt(math.log(1 / eps)))}
        print(f"  fit: d = {a:.3f} alpha/Delta + {b:.0f};  slope/sqrt(ln 1/eps) = {a/math.sqrt(math.log(1/eps)):.3f}"
              f"  (model constant {4.7571/math.sqrt(math.log(100)):.3f})", flush=True)
    c0, e0, ecore = identity_and_e0()
    out["bodipy_core_cas12_11"] = {"c0_Ha": c0, "E0_Ha": e0, "ecore_Ha": ecore, "abs_E0_minus_c0_Ha": abs(e0 - c0)}
    print(out["bodipy_core_cas12_11"], flush=True)
    path = os.path.join(ROOT, "paper", "quantum", "results", "degree_check.json")
    json.dump(out, open(path, "w"), indent=1)


if __name__ == "__main__":
    main()
