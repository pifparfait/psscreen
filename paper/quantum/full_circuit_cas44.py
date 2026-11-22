"""End-to-end simulation of one call of the decision circuit in CAS(4,4): trial state, S1 filter, H_SO block
encoding with normalization alpha_SO, and the golden-rule QSVT filter.

decision_cost.py assumes that one call succeeds with probability p = gamma^2 k_ISC / K', K' = (2 pi/hbar)
alpha_SO^2 rho_max, where gamma^2 is the success probability of the S1 preparation. ae_demo.py and isc_filtered.py
check the filter alone, starting from chi = H_SO|S1>/||chi|| prepared exactly. Here the whole call is one circuit on
the 8 spin-orbital qubits of a CAS(4,4) (wire 2p = p alpha, 2p+1 = p beta) plus flags:
  1. trial state: the spin-adapted HOMO->LUMO singlet configuration (no triplet component by construction);
  2. S1 eigenstate filter: QSVT of an even polynomial ~ amp_S exp(-(E - E_S1)^2 / (2 sigma^2)) on a block encoding
     of A = (H - E_c)/alpha (one block-encoding qubit, one qubit for the real part);
  3. H_SO / alpha_SO as a linear combination of unitaries: after the orbital rotation that diagonalizes the
     spin-orbital matrix M, H_SO = -(1/2) sum_k mu_k Z_k, so PREP loads sqrt(|mu_k|/sum|mu|) on a 3-qubit index
     register and SELECT applies -sign(mu_k) Z_k; the rotation is the Fock-space unitary exp(sum_pq K_pq a+_p a_q)
     with K = log u, u the eigenvectors of M;
  4. golden-rule filter: QSVT of ~ 0.98 sqrt(rho_MLJ/rho_max) on the same block encoding (two more flag qubits).
Success = all flags in |0>. The prediction is p = |<S1|phi>|^2 g_S(E_S1)^2 0.98^2 k_ISC/K' (up to polynomial errors
and the leakage of the S1 filter), where k_ISC is the exact golden rule of the same Hamiltonian. As in
isc_filtered.py, the block encodings of H are explicit dilations (built from the matrix; on hardware they would be
the THC walk), and the number-conserving encoding is simulated in the 4-electron sector with the other sectors
padded.

Usage (repository root): python paper/quantum/full_circuit_cas44.py [NAME]
Writes paper/quantum/results/NAME.full_circuit_cas44.json.
"""
import json
import math
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "paper", "quantum"))
import numpy as np  # noqa: E402
import pennylane as qml  # noqa: E402
from numpy.polynomial import chebyshev as C  # noqa: E402
from pyscf import dft, gto  # noqa: E402
from scipy.linalg import expm, logm  # noqa: E402

from alpha_so import spin_orbital_matrix  # noqa: E402
from psscreen.backends import _ecp  # noqa: E402
from psscreen.geometry import read_xyz, to_pyscf_atom  # noqa: E402
from psscreen.isc_filtered import filter_poly  # noqa: E402
from psscreen.qsp_angles import qsp_phases_fpi, to_pennylane_qsvt  # noqa: E402
from psscreen.quantum_pl import HARTREE_EV, active_space, eigenstates, qubit_hamiltonian, qubit_onebody  # noqa: E402
from psscreen.rates import CM1_TO_EV, HBAR_EVS, MLJParams, fcwd_mlj  # noqa: E402

NCAS, NQ = 4, 8
HA_CM = 219474.6313705
P = MLJParams()


def fermi_matrix(K):
    """Dense matrix of sum_pq K_pq a+_p a_q on the 8 qubits (Jordan-Wigner, PennyLane wire order)."""
    fs = qml.fermi.FermiSentence({})
    for p in range(NQ):
        for q in range(NQ):
            if abs(K[p, q]) > 1e-14:
                fs += complex(K[p, q]) * qml.FermiWord({(0, p): "+", (1, q): "-"})
    return qml.matrix(qml.jordan_wigner(fs), wire_order=range(NQ))


def number_mask():
    return np.array([bin(i).count("1") for i in range(2 ** NQ)])


def dilation(Hfull, e_c, alpha, nmask, ne):
    """512x512 unitary with block A = (H - e_c)/alpha on the ne-electron sector, 0.999 elsewhere."""
    dim = 2 ** NQ
    A = np.eye(dim, dtype=complex) * 0.999
    sel = np.where(nmask == ne)[0]
    A[np.ix_(sel, sel)] = (Hfull[np.ix_(sel, sel)] - e_c * np.eye(len(sel))) / alpha
    A = 0.5 * (A + A.conj().T)
    ea, va = np.linalg.eigh(A)
    assert np.abs(ea).max() < 1
    Sq = (va * np.sqrt(1 - ea ** 2)) @ va.conj().T
    return np.block([[A, Sq], [Sq, -A]])


def even_poly(f, tol, max_degree=6000):
    deg = 40
    xs = np.linspace(-1, 1, 40001)
    fx = f(xs)
    while True:
        deg += deg % 2
        c = C.chebinterpolate(f, deg)
        err = float(np.max(np.abs(C.chebval(xs, c) - fx)))
        if err < tol or deg >= max_degree:
            return c[0::2], deg, err
        deg = int(deg * 1.5)


def main(name):
    t0 = time.time()
    path = os.path.join(ROOT, "validation", "runs", f"{name}.xyz")
    sym, xyz = read_xyz(path)
    atom = to_pyscf_atom(sym, xyz)
    mol = gto.M(atom=atom, basis="def2-svp", ecp=_ecp(atom, "def2-svp"), verbose=0)
    mf = dft.RKS(mol, xc="b3lyp").density_fit()
    mf.grids.level = 2
    mf.kernel()
    asp = active_space(mf, NCAS, NCAS, "def2-svp")
    ne = asp["nelec"]
    Hq = qubit_hamiltonian(asp)
    st = eigenstates(asp, Hq)
    Hfull = st["Hs"].toarray()
    nmask = number_mask()
    E0 = st["E0"]
    E_S1 = E0 + st["E_S"][1] / HARTREE_EV
    S1 = st["S1"]

    # spin-orbital H_SO and its diagonal form
    Hso = qml.matrix(qml.simplify(qubit_onebody(soc=asp["soc"], n=NCAS)), wire_order=range(NQ))
    herm_err = float(np.abs(Hso - Hso.conj().T).max())
    Hso = 0.5 * (Hso + Hso.conj().T)     # the integrals are antisymmetric only to ~1e-3 relative; keep the Hermitian part
    M = spin_orbital_matrix(asp["soc"])
    M = 0.5 * (M + M.conj().T)
    so_err = float(np.abs(fermi_matrix(M) - Hso).max())
    print(f"|M (Fock) - qubit H_SO| = {so_err:.2e} Ha (max element {np.abs(Hso).max():.2e})", flush=True)
    # qml.simplify drops Pauli terms below ~1e-8, hence a difference of that size; the exact matrix is used below
    assert so_err < 1e-7, "spin-orbital matrix and qubit H_SO disagree"
    Hso = fermi_matrix(M)
    mu, u = np.linalg.eigh(M)
    alpha_so = 0.5 * float(np.abs(mu).sum())
    UF = expm(fermi_matrix(logm(u)))
    Hdiag = sum(mu[k] * fermi_matrix(np.outer(np.eye(NQ)[k], np.eye(NQ)[k])) for k in range(NQ))
    rot_err = float(np.abs(UF @ Hdiag @ UF.conj().T - Hso).max())
    w = np.abs(mu) / np.abs(mu).sum()
    sgn = -np.sign(mu)

    # trial state: (E_{LUMO,HOMO}^alpha + E^beta)|HF>/sqrt(2), HOMO = orbital 1, LUMO = orbital 2
    hf = np.zeros(2 ** NQ, dtype=complex)
    hf[int("11110000", 2)] = 1.0
    exc = fermi_matrix(np.eye(NQ)[:, [4]] @ np.eye(NQ)[[2], :] + np.eye(NQ)[:, [5]] @ np.eye(NQ)[[3], :])
    phi = exc @ hf
    phi /= np.linalg.norm(phi)
    gamma2 = float(abs(np.vdot(S1, phi)) ** 2)

    # block encoding of H and the two polynomials
    sel = np.where(nmask == ne)[0]
    ev = np.linalg.eigvalsh(Hfull[np.ix_(sel, sel)])
    e_c = ev[0] - 1e-3
    alpha = 1.01 * (ev[-1] - e_c)
    UA = dilation(Hfull, e_c, alpha, nmask, ne)
    Es = st["E0"] + np.asarray(st["E_S"]) / HARTREE_EV
    gap = min(abs(Es[0] - E_S1), abs(Es[2] - E_S1))
    sigma = gap / 3
    amp_s = 0.98

    def g(x):
        E = alpha * np.abs(np.asarray(x, float)) + e_c
        return amp_s * np.exp(-0.5 * ((E - E_S1) / sigma) ** 2)

    cS, degS, errS = even_poly(g, 1e-4)
    angS = to_pennylane_qsvt(qsp_phases_fpi(cS, degS)[0])
    f, cI, degI, rho_max, errI = filter_poly(alpha, e_c, E_S1, P, tol=2e-3)
    angI = to_pennylane_qsvt(qsp_phases_fpi(cI, degI)[0])
    print(f"gamma2 {gamma2:.4f}  S1 filter degree {degS} (err {errS:.1e}), ISC filter degree {degI} "
          f"(err {errI:.1e}), alpha_SO {alpha_so*HA_CM:.1f} cm-1, rotation err {rot_err:.1e}", flush=True)

    # wires: 0 c1, 1 a1, 2 c2, 3 a2, 4-6 index, 7-14 system
    c1, a1, c2, a2 = 0, 1, 2, 3
    idx = [4, 5, 6]
    sysw = list(range(7, 15))
    prep = np.sqrt(w)
    selects = [np.diag([sgn[k] * (1 if not (i >> (NQ - 1 - k)) & 1 else -1) for i in range(2 ** NQ)])
               for k in range(NQ)]
    dev = qml.device("default.qubit", wires=15)

    def qsvt_re(ctl, anc, angles):
        qml.Hadamard(ctl)
        for s, cv in ((1, 0), (-1, 1)):
            proj = [qml.PCPhase(s * a_, dim=2 ** NQ, wires=[anc] + sysw) for a_ in angles]
            qml.ctrl(qml.QSVT(qml.QubitUnitary(UA, wires=[anc] + sysw), proj), control=ctl, control_values=cv)
        qml.Hadamard(ctl)

    @qml.qnode(dev)
    def call(stage):
        qml.StatePrep(phi, wires=sysw)
        qsvt_re(c1, a1, angS)
        if stage >= 2:
            qml.QubitUnitary(UF.conj().T, wires=sysw)
            qml.MottonenStatePreparation(prep, wires=idx)
            for k in range(NQ):
                qml.ctrl(qml.QubitUnitary(selects[k], wires=sysw), control=idx,
                         control_values=[int(b) for b in format(k, "03b")])
            qml.adjoint(qml.MottonenStatePreparation)(prep, wires=idx)
            qml.QubitUnitary(UF, wires=sysw)
        if stage >= 3:
            qsvt_re(c2, a2, angI)
        return qml.state()

    def good_prob(psi, flags):
        psi = psi.reshape([2] * 15)
        sl = tuple(0 if q in flags else slice(None) for q in range(15))
        return float(np.sum(np.abs(psi[sl]) ** 2)), psi[sl]

    p1, v1 = good_prob(call(1), [c1, a1, c2, a2] + idx)
    s1_fid = float(abs(np.vdot(S1, v1.reshape(-1))) ** 2 / p1)
    p2, _ = good_prob(call(2), [c1, a1, c2, a2] + idx)
    t1 = time.time()
    p3, _ = good_prob(call(3), [c1, a1, c2, a2] + idx)
    t_sim = time.time() - t1

    # exact golden rule on the same Hamiltonian, all M_S sectors
    chi = Hso @ S1
    e_all, v_all = np.linalg.eigh(Hfull[np.ix_(sel, sel)])
    ov = np.abs(v_all.conj().T @ chi[sel]) ** 2
    dE = (e_all - E_S1) * HARTREE_EV
    rho = np.array([fcwd_mlj(float(x), P) for x in dE])
    pref = 2 * math.pi / HBAR_EVS * HARTREE_EV ** 2
    k_exact = pref * float(np.sum(ov * rho))
    Kp = 2 * math.pi / HBAR_EVS * (alpha_so * HA_CM * CM1_TO_EV) ** 2 * rho_max
    gS = float(g((E_S1 - e_c) / alpha))
    p_pred = gamma2 * gS ** 2 * 0.98 ** 2 * k_exact / Kp
    out = {"name": name, "active_space": "CAS(4,4), B3LYP/def2-SVP orbitals", "gamma2": gamma2,
           "S1_filter": {"degree": degS, "poly_err": errS, "sigma_eV": sigma * HARTREE_EV, "gap_eV": gap * HARTREE_EV,
                         "success": p1, "fidelity_with_S1": s1_fid},
           "alpha_so_cm": alpha_so * HA_CM, "rotation_err_Ha": rot_err, "so_matrix_err_Ha": so_err, "hso_antihermitian_part_Ha": herm_err,
           "p_after_HSO": p2, "p_after_HSO_pred": p1 * s1_fid * float(np.vdot(chi, chi).real) / alpha_so ** 2,
           "isc_filter": {"degree": degI, "poly_err": errI, "alpha_Ha": alpha},
           "p_call": p3, "p_pred": p_pred, "ratio_call_to_pred": p3 / p_pred, "k_exact": k_exact, "K_prime": Kp,
           "k_from_call": p3 / (gamma2 * gS ** 2 * 0.98 ** 2) * Kp, "seconds_sim": t_sim,
           "seconds": time.time() - t0}
    print(json.dumps(out, indent=1), flush=True)
    json.dump(out, open(os.path.join(ROOT, "paper", "quantum", "results", f"{name}.full_circuit_cas44.json"), "w"),
              indent=1)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "BODIPY_core")
