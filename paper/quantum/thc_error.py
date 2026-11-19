"""Error of the tensor-hypercontraction (THC) factorization on the excitation energies that the screening needs.

The resource estimates of the paper assume a THC rank of 22 (2N) for the (12e,11o) BODIPY-core Hamiltonian and an
assumed one-norm. Here the two-electron integrals of that Hamiltonian (canonical B3LYP window, def2-SVP, the same
orbitals as paper/quantum/resource_breakdown.py) are fitted to the THC form
    (pq|rs) ~ sum_{mu,nu} X_{p mu} X_{q mu} Z_{mu nu} X_{r nu} X_{s nu}
by least squares (L-BFGS with analytic gradients, several random starts, best kept), for ranks 2N, 3N and 4N. For each
rank the script reports the integral error, the THC one-norm lambda = lambda_T + lambda_Z in the convention of Lee et al.
(PRX Quantum 2, 030305 (2021)), and the CASCI energies of S1, T1 and T2 with the exact and the THC integrals, i.e. the
error that THC adds to the S1-T1 and S1-T2 gaps compared with the accuracy budget.
Usage (repository root): python paper/quantum/thc_error.py [NAME] [--ranks 22 33 44] [--starts 6]
Writes paper/quantum/results/NAME.thc_error.json.
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
ap.add_argument("--ranks", type=int, nargs="+", default=[22, 33, 44])
ap.add_argument("--starts", type=int, default=6)
ap.add_argument("--maxiter", type=int, default=20000)
a = ap.parse_args()
OUT = os.path.join(ROOT, "paper", "quantum", "results", f"{a.name}.thc_error.json")

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


def loss_grad(x, M):
    X, Z = unpack(x, M)
    T = np.einsum("pm,qm->mpq", X, X).reshape(M, N * N)
    R = T.T @ Z @ T - Vm
    L = (R ** 2).sum()
    GT = 2 * (Z @ T @ R.T + Z.T @ T @ R)
    GZ = 2 * T @ R @ T.T
    GT = GT.reshape(M, N, N)
    GX = np.einsum("mpq,qm->pm", GT, X) + np.einsum("mqp,qm->pm", GT, X)
    GZu = (GZ + GZ.T - np.diag(np.diag(GZ)))[np.triu_indices(M)]
    return L, np.concatenate([GX.ravel(), GZu])


def fit(M, seed):
    rng = np.random.default_rng(seed)
    # start: X from the leading eigenvectors of the pair-density matrix, then noise
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
    r = minimize(loss_grad, x0, args=(M,), jac=True, method="L-BFGS-B",
                 options={"maxiter": a.maxiter, "maxfun": 4 * a.maxiter, "ftol": 1e-16, "gtol": 1e-12})
    return r


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
       "exact": states(V), "ranks": {}}
print("exact", res["exact"], flush=True)
for M in a.ranks:
    best = None
    for k in range(a.starts):
        r = fit(M, k)
        if best is None or r.fun < best.fun:
            best = r
        print(f"rank {M} start {k}: loss {r.fun:.3e} ({r.nit} it)", flush=True)
    X, Z = unpack(best.x, M)
    Vt = thc_eri(X, Z)
    st = states(Vt)
    lt, lz = one_norm(h1, V, X, Z)
    res["ranks"][str(M)] = {
        "eri_max_abs_err_Ha": float(abs(Vt - V).max()), "eri_rms_err_Ha": float(np.sqrt(((Vt - V) ** 2).mean())),
        "states": st, "lambda_T_Ha": lt, "lambda_Z_Ha": lz, "lambda_Ha": lt + lz,
        "gap_errors_eV": {k: st[k] - res["exact"][k] for k in st}}
    print(M, json.dumps(res["ranks"][str(M)]), flush=True)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(res, open(OUT, "w"), indent=1)
res["time_s"] = time.time() - t0
json.dump(res, open(OUT, "w"), indent=1)
print("done", round(time.time() - t0), "s")
