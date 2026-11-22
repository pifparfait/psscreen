"""Success probability p of the gap-resolved filter from one and the same CAS(12,11) Hamiltonian.

decision_cost.py combines the pipeline rate k_ISC (TDDFT/DeltaSCF) with the block-encoding norm alpha_SO of the
CAS(12,11) spin-orbit operator, p = gamma^2 k_ISC / K'. Here both come from the active space: exact CASCI singlets
and triplets of the (12e,11o) canonical window with B3LYP/def2-SVP orbitals, couplings from the same calibrated
spin-orbit integrals, and the same MLJ line shape:
    V_n^2 = sum_M |<S1|H_SO|T_n,M>|^2 = sum_k |sum_pq h^k_pq rho^n_pq|^2,
    rho^n_pq = <S1| (a+_pa a_qa - a+_pb a_qb)/2 |T_n, M=0>   (Wigner-Eckart: the sum over M equals the sum over k;
    T_n is solved at M = 1 and lowered with S-),
    k_CAS = sum_n (2 pi/hbar) V_n^2 rho_MLJ(E_Tn - E_S1),   p_CAS = gamma^2 k_CAS / K'.
This is the success probability the filter would actually see on that Hamiltonian, to compare with the
pipeline-target p used for the costs.

Usage (repository root): python paper/quantum/cas_consistent_p.py NAME [NAME ...] [--ntrip 6]
Writes paper/quantum/results/NAME.cas_p.json.
"""
import argparse
import json
import math
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "paper", "quantum"))
import numpy as np  # noqa: E402
from pyscf import dft, fci, gto, lib  # noqa: E402

from psscreen.backends import _ecp  # noqa: E402
from psscreen.geometry import read_xyz, to_pyscf_atom  # noqa: E402
from psscreen.quantum_pl import active_space  # noqa: E402
from psscreen.rates import CM1_TO_EV, HBAR_EVS, MLJParams, fcwd_mlj  # noqa: E402
from alpha_so import XYZ, spin_orbital_matrix  # noqa: E402

EV = 27.211386245988
HA_CM = 219474.6313705
GAMMA2 = 0.5
P = MLJParams()
RHO_MAX = max(fcwd_mlj(float(x), P) for x in np.linspace(-1.5, 0.5, 2001))


def solve(h1, h2, N, nelec, ss, nroots, ecore):
    s = fci.direct_spin1.FCI()
    s.nroots, s.conv_tol = nroots, 1e-10
    s = fci.addons.fix_spin_(s, ss=ss)
    e, c = s.kernel(h1, h2, N, nelec, ecore=ecore, max_cycle=400)
    return np.asarray(e), c


def lower(ci, N):
    """S- = sum_p a+_{p beta} a_{p alpha} applied to an M = 1 vector of (7, 5) electrons, normalized."""
    out = 0
    for p in range(N):
        t = fci.addons.des_a(ci, N, (7, 5), p)
        out = out + fci.addons.cre_b(t, N, (6, 5), p)
    return out / np.linalg.norm(out)


def run(name, ntrip):
    t0 = time.time()
    sym, xyz = read_xyz(os.path.join(ROOT, XYZ[name]))
    atom = to_pyscf_atom(sym, xyz)
    mol = gto.M(atom=atom, basis="def2-svp", ecp=_ecp(atom, "def2-svp"), verbose=0)
    mf = dft.RKS(mol, xc="b3lyp").density_fit()
    mf.grids.level = 2
    mf.kernel()
    asp = active_space(mf, 11, 12, "def2-svp")
    h1, h2, N, ec = asp["h1"], asp["h2"], asp["ncas"], asp["ecore"]
    soc = asp["soc"]                                   # complex, H_SO = sum_k h^k_pq s_k (spin-orbital form)
    mu = np.linalg.eigvalsh(0.5 * (spin_orbital_matrix(soc) + spin_orbital_matrix(soc).conj().T))
    alpha = 0.5 * float(np.abs(mu).sum())
    eS, cS = solve(h1, h2, N, (6, 6), 0, 3, ec)
    eT, cT1 = solve(h1, h2, N, (7, 5), 2, ntrip, ec)
    cT = [lower(c, N) for c in cT1]              # M = 0 components, S-|1,1> = sqrt(2)|1,0>
    s1 = cS[1]
    chans = []
    for n in range(ntrip):
        da, db = fci.direct_spin1.trans_rdm1s(s1, cT[n], N, (6, 6))
        rho = 0.5 * (da - db)
        v2 = sum(abs(np.einsum("pq,pq->", soc[k], rho.T)) ** 2 for k in range(3))
        V = math.sqrt(v2) * HA_CM
        dE = (eT[n] - eS[1]) * EV
        k = 2 * math.pi / HBAR_EVS * (V * CM1_TO_EV) ** 2 * fcwd_mlj(dE, P)
        chans.append({"T": n + 1, "E_T_eV": float((eT[n] - eS[0]) * EV), "dE_eV": float(dE), "V_cm": V, "k": k})
    kcas = sum(c["k"] for c in chans)
    Kp = 2 * math.pi / HBAR_EVS * (alpha * HA_CM * CM1_TO_EV) ** 2 * RHO_MAX
    out = {"name": name, "active_space": "(12e,11o) canonical window, B3LYP/def2-SVP orbitals, exact CASCI",
           "E_S1_eV": float((eS[1] - eS[0]) * EV), "E_S2_eV": float((eS[2] - eS[0]) * EV),
           "S1_T1_eV": float((eS[1] - eT[0]) * EV), "alpha_so_cm": alpha * HA_CM, "K_prime": Kp,
           "channels": chans, "k_isc_cas": kcas, "p_cas": GAMMA2 * kcas / Kp,
           "chi_norm_cm_triplets_computed": math.sqrt(sum(c["V_cm"] ** 2 for c in chans)),
           "seconds": time.time() - t0}
    path = os.path.join(ROOT, "paper", "quantum", "results", f"{name}.cas_p.json")
    json.dump(out, open(path, "w"), indent=1)
    print(name, json.dumps({k: v for k, v in out.items() if k != "channels"}), flush=True)
    for c in chans:
        print("  ", c, flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("names", nargs="+")
    ap.add_argument("--ntrip", type=int, default=6)
    a = ap.parse_args()
    lib.num_threads(4)
    for n in a.names:
        run(n, a.ntrip)
