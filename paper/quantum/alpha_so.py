"""Block-encoding normalization alpha_SO of the one-body spin-orbit operator in the CAS(12,11) active space.

H_SO = sum_{PQ} M_PQ a_P^dagger a_Q over the 2N active spin-orbitals (M Hermitian, from the calibrated one-electron
SOC integrals of psscreen.soc). After one orbital rotation that diagonalizes M, H_SO = sum_k mu_k n_k
= (1/2) tr M - (1/2) sum_k mu_k Z_k, so the LCU normalization of the traceless part is
    alpha_SO = (1/2) sum_k |mu_k|.
The success probability of one call of the gap-resolved filter is then p_tot = <chi|P^2|chi> / alpha_SO^2, with
chi = H_SO|S1>, and k_ISC = (2 pi/hbar) alpha_SO^2 rho_max p_tot.

Usage (from the repository root): python paper/quantum/alpha_so.py [NAME ...]
Writes paper/quantum/results/alpha_so.json (merged with existing entries).
"""
import json
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
import numpy as np  # noqa: E402
from pyscf import dft, gto, lib  # noqa: E402

from psscreen.backends import _ecp  # noqa: E402
from psscreen.geometry import read_xyz, to_pyscf_atom  # noqa: E402
from psscreen.quantum_pl import active_space  # noqa: E402

HA_CM = 219474.6313705
NCAS, NELE = 11, 12
XYZ = {
    **{n: f"validation/runs/{n}.xyz" for n in ["BODIPY_core", "BDP", "I-BDP", "I2-BDP", "I3-BDP", "azaBDP_Ph4"]},
    **{n: f"validation/runs/set3/{n}_c0.xyz" for n in ["BDP-0", "I-BDP-0", "I2-BDP-0"]},
    **{n: f"screening/runs/candidates_round1/{n}.xyz"
       for n in ["Cl2-BDP", "Br2-BDP", "I1-BDP", "azaBDP_Br2", "azaBDP_BrPh2", "azaBDP_IPh2"]},
}
PAULI = [np.array([[0, 1], [1, 0]]) / 2, np.array([[0, -1j], [1j, 0]]) / 2, np.array([[1, 0], [0, -1]]) / 2]


def spin_orbital_matrix(soc):
    """M[(p,s),(q,t)] = sum_k h^k_pq (sigma_k/2)_st, the same convention as quantum_pl.qubit_onebody."""
    n = soc.shape[1]
    hs = sum(np.einsum("pq,st->spqt", soc[k], PAULI[k]) for k in range(3))   # [s,p,q,t]
    return np.einsum("spqt->psqt", hs).reshape(2 * n, 2 * n)


def main(names):
    out_path = os.path.join(ROOT, "paper", "quantum", "results", "alpha_so.json")
    res = json.load(open(out_path)) if os.path.exists(out_path) else {}
    for name in names:
        t0 = time.time()
        sym, xyz = read_xyz(os.path.join(ROOT, XYZ[name]))
        atom = to_pyscf_atom(sym, xyz)
        mol = gto.M(atom=atom, basis="def2-svp", ecp=_ecp(atom, "def2-svp"), verbose=0)
        mf = dft.RKS(mol, xc="b3lyp").density_fit()
        mf.grids.level = 2
        mf.kernel()
        asp = active_space(mf, NCAS, NELE, "def2-svp")
        M = spin_orbital_matrix(asp["soc"])
        herm = float(np.abs(M - M.conj().T).max())
        mu = np.linalg.eigvalsh(0.5 * (M + M.conj().T))
        alpha = 0.5 * float(np.sum(np.abs(mu)))
        res[name] = {"alpha_so_cm": alpha * HA_CM, "max_abs_mu_cm": float(np.abs(mu).max()) * HA_CM,
                     "trace_cm": float(np.real(np.trace(M))) * HA_CM, "hermiticity_err": herm,
                     "converged": bool(mf.converged), "natm": mol.natm, "seconds": time.time() - t0}
        print(name, json.dumps(res[name]), flush=True)
        json.dump(res, open(out_path, "w"), indent=1)


if __name__ == "__main__":
    lib.num_threads(4)
    main(sys.argv[1:] or list(XYZ))
