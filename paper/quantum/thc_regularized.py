"""THC fit with a one-norm penalty: lambda and the gap error of the regularized rank-33 factorization.

Same Hamiltonian, fit and CASCI states as paper/quantum/thc_error.py; the least-squares loss gets the penalty
rho * lambda_Z^2 of Lee et al. (PRX Quantum 2, 030305 (2021)), and each rho is warm-started from the previous one,
starting from the best unregularized fit. Usage (repository root):
    python paper/quantum/thc_regularized.py [NAME] [--ranks 33] [--starts 4] [--rhos 0 1e-6 1e-5 1e-4 1e-3]
Writes paper/quantum/results/NAME.thc_regularized.json.
"""
import argparse
import json
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
import numpy as np  # noqa: E402
from pyscf import dft, fci, gto  # noqa: E402
from scipy.optimize import minimize  # noqa: E402

from psscreen.backends import _ecp  # noqa: E402
from psscreen.geometry import read_xyz, to_pyscf_atom  # noqa: E402
from psscreen.quantum_pl import active_space  # noqa: E402

EV = 27.211386245988
ap = argparse.ArgumentParser()
ap.add_argument("name", nargs="?", default="BODIPY_core")
ap.add_argument("--ranks", type=int, nargs="+", default=[33])
ap.add_argument("--starts", type=int, default=4)
ap.add_argument("--maxiter", type=int, default=20000)
ap.add_argument("--rhos", type=float, nargs="+", default=[0.0, 1e-6, 1e-5, 1e-4, 1e-3])
a = ap.parse_args()
OUT = os.path.join(ROOT, "paper", "quantum", "results", f"{a.name}.thc_regularized.json")
EPS = 1e-8

sym, xyz = read_xyz(os.path.join(ROOT, "validation", "runs", f"{a.name}.xyz"))
atom = to_pyscf_atom(sym, xyz)
mol = gto.M(atom=atom, basis="def2-svp", ecp=_ecp(atom, "def2-svp"), verbose=0)
mf = dft.RKS(mol, xc="b3lyp").density_fit()
mf.grids.level = 2
mf.kernel()
asp = active_space(mf, 11, 12, "def2-svp")
h1, V, N = asp["h1"], asp["h2"], asp["ncas"]          # V[p,q,r,s] = (pq|rs)
Vm = V.reshape(N * N, N * N)


def unpack(x, M):
    X = x[:N * M].reshape(N, M)
    Zu = x[N * M:]
    Z = np.zeros((M, M))
    Z[np.triu_indices(M)] = Zu
    return X, Z + np.triu(Z, 1).T


def loss_grad(x, M, rho=0.0):
    X, Z = unpack(x, M)
    T = np.einsum("pm,qm->mpq", X, X).reshape(M, N * N)
    R = T.T @ Z @ T - Vm
    L = (R ** 2).sum()
    GT = 2 * (Z @ T @ R.T + Z.T @ T @ R)
    GZ = 2 * T @ R @ T.T
    GT = GT.reshape(M, N, N)
    GX = np.einsum("mpq,qm->pm", GT, X) + np.einsum("mqp,qm->pm", GT, X)
    if rho:
        # one-norm penalty rho * lambda_Z^2 (Lee et al. 2021), |Z| smoothed as sqrt(Z^2 + EPS^2)
        nrm = (X ** 2).sum(0)
        A = np.sqrt(Z ** 2 + EPS ** 2)
        lz = 0.5 * (A * np.outer(nrm, nrm)).sum()
        L += rho * lz ** 2
        GZ = GZ + rho * 2 * lz * 0.5 * (Z / A) * np.outer(nrm, nrm)
        GX = GX + rho * 2 * lz * 2 * X * (A @ nrm)[None, :]
    GZu = (GZ + GZ.T - np.diag(np.diag(GZ)))[np.triu_indices(M)]
    return L, np.concatenate([GX.ravel(), GZu])


def fit(M, seed, rho=0.0, x0=None):
    if x0 is None:
        rng = np.random.default_rng(seed)
        w, U = np.linalg.eigh(Vm)
        X0 = np.zeros((N, M))
        for m in range(M):
            u = U[:, -1 - (m % (N * N))].reshape(N, N)
            e, v = np.linalg.eigh((u + u.T) / 2)
            X0[:, m] = v[:, np.argmax(abs(e))] * np.sqrt(abs(e).max())
        X0 += 0.1 * rng.standard_normal(X0.shape)
        T = np.einsum("pm,qm->mpq", X0, X0).reshape(M, N * N)
        Z0 = np.linalg.lstsq(T @ T.T, T @ Vm @ T.T @ np.linalg.pinv(T @ T.T), rcond=None)[0]
        Z0 = (Z0 + Z0.T) / 2
        x0 = np.concatenate([X0.ravel(), Z0[np.triu_indices(M)]])
    return minimize(loss_grad, x0, args=(M, rho), jac=True, method="L-BFGS-B",
                    options={"maxiter": a.maxiter, "maxfun": 4 * a.maxiter, "ftol": 1e-16, "gtol": 1e-12})


def thc_eri(X, Z):
    return np.einsum("pm,qm,mn,rn,sn->pqrs", X, X, Z, X, X, optimize=True)


def one_norm(h1, V, X, Z):
    """lambda_T + lambda_Z (Lee et al. 2021, Eqs. 12-13 conventions, normalized THC leaf vectors)."""
    T = h1 - 0.5 * np.einsum("prrq->pq", V) + np.einsum("pqrr->pq", V)
    lam_T = np.abs(np.linalg.eigvalsh(T)).sum()
    nrm = (X ** 2).sum(0)
    lam_Z = 0.5 * (np.abs(Z) * np.outer(nrm, nrm)).sum()
    return float(lam_T), float(lam_Z)


def states(h2):
    out = {}
    for spin, n in ((0, 3), (2, 3)):
        s = fci.direct_spin1.FCI()
        s.spin, s.nroots = spin, n
        s = fci.addons.fix_spin_(s, ss=spin / 2 * (spin / 2 + 1))
        e, _ = s.kernel(h1, h2, N, ((12 + spin) // 2, (12 - spin) // 2), ecore=asp["ecore"], max_cycle=400)
        out["S" if spin == 0 else "T"] = np.asarray(e)
    S, T = out["S"], out["T"]
    return {"S1_eV": float((S[1] - S[0]) * EV), "T1_eV": float((T[0] - S[0]) * EV), "T2_eV": float((T[1] - S[0]) * EV),
            "S1_T1_eV": float((S[1] - T[0]) * EV), "T2_S1_eV": float((T[1] - S[1]) * EV)}


t0 = time.time()
res = {"name": a.name, "active_space": "(12e,11o) canonical window, B3LYP/def2-SVP orbitals", "N": N,
       "penalty": "rho * lambda_Z^2, warm-started from the unregularized fit", "exact": states(V), "ranks": {}}
print("exact", res["exact"], flush=True)
for M in a.ranks:
    best = None
    for k in range(a.starts):
        r = fit(M, k)
        if best is None or r.fun < best.fun:
            best = r
        print(f"rank {M} start {k}: loss {r.fun:.3e} ({r.nit} it)", flush=True)
    res["ranks"][str(M)] = {}
    x = best.x
    for rho in a.rhos:
        r = best if rho == 0 else fit(M, 0, rho, x)
        x = r.x
        X, Z = unpack(r.x, M)
        Vt = thc_eri(X, Z)
        st = states(Vt)
        lt, lz = one_norm(h1, V, X, Z)
        d = {"eri_max_abs_err_Ha": float(abs(Vt - V).max()), "eri_rms_err_Ha": float(np.sqrt(((Vt - V) ** 2).mean())),
             "states": st, "lambda_T_Ha": lt, "lambda_Z_Ha": lz, "lambda_Ha": lt + lz,
             "gap_errors_eV": {k: st[k] - res["exact"][k] for k in st}}
        res["ranks"][str(M)][f"{rho:g}"] = d
        print(M, rho, f"lambda {lt + lz:.2f}", {k: round(v, 4) for k, v in d["gap_errors_eV"].items()}, flush=True)
        os.makedirs(os.path.dirname(OUT), exist_ok=True)
        json.dump(res, open(OUT, "w"), indent=1)
res["time_s"] = time.time() - t0
json.dump(res, open(OUT, "w"), indent=1)
print("done", round(time.time() - t0), "s")
