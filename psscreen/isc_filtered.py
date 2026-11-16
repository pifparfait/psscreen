"""Gap-resolved ISC observable: Fermi's golden rule as a QSVT spectral filter (PennyLane, simulator).

The evolution proxy of arXiv:2512.15889, |<T1,M| exp(-i H_SO t) |S1>|^2, estimates the S1-T1 spin-orbit
coupling but contains no energy gap and no Franck-Condon factor, and it looks at one triplet only.
The golden-rule rate over all triplets can be written as one expectation value,

    k_ISC = (2 pi / hbar) <chi| rho(H - E_S1) |chi>,      |chi> = H_SO |S1>,

where rho is the Franck-Condon weighted density of states (here the Marcus-Levich-Jortner form of
psscreen.rates) evaluated on the spectrum of the electronic Hamiltonian H. Because H_SO couples S1
only to triplets, the sum over eigenstates of H in |chi> runs over every triplet T_n, M with its own
gap. This is the same structure as the absorption-window observable of the paper, with the dipole
operator replaced by H_SO and the window indicator replaced by rho, so it can be computed with the same
machinery: block-encode A = (H - E_c)/alpha, apply a QSVT polynomial P(A) ~ sqrt(rho/rho_max), and read
the success probability p = ||P(A) chi||^2 / ||chi||^2. Then

    k_ISC = (2 pi / hbar) ||chi||^2 rho_max p.

The electronic part (energies, couplings, all triplets) is quantum; the vibronic line shape rho is a
classical model with assumed parameters, as in the classical pipeline.

As in quantum_pl.threshold_projection, the simulation builds an explicit block encoding in a generic
(randomly rotated) basis of each M_S sector, so no eigenbasis information enters the circuit. S1 is
prepared exactly; E_S1 is an input (from phase estimation or the threshold projection).
"""
from __future__ import annotations

import numpy as np
import pennylane as qml

from .qsp_angles import qsp_phases_fpi, to_pennylane_qsvt
from .quantum_pl import HARTREE_EV, _sector, qubit_onebody
from .rates import CM1_TO_EV, HBAR_EVS, MLJParams, fcwd_mlj
from .soc import CM1_PER_HARTREE


def _rho_vec(dE_eV, p):
    return np.array([fcwd_mlj(float(x), p) for x in np.atleast_1d(dE_eV)])


def filter_poly(alpha, e_c, e_s1, p, amp=0.98, tol=2e-3, max_degree=8000):
    """Even Chebyshev polynomial ~ amp*sqrt(rho(alpha*|x| + e_c - e_s1)/rho_max) on [-1, 1] (energies in Ha).
    Returns (f, even Chebyshev coefficients, degree, rho_max, max approximation error)."""
    from numpy.polynomial import chebyshev as C
    grid = np.linspace(-1.5, 0.5, 2001)
    rho_max = float(np.max(_rho_vec(grid, p)))

    def f(x):
        dE = (alpha * np.abs(np.asarray(x, dtype=float)) + e_c - e_s1) * HARTREE_EV
        return amp * np.sqrt(np.clip(_rho_vec(dE.ravel(), p).reshape(np.shape(dE)), 0, None) / rho_max)

    # the filter width in x sets the degree; start from ~8/width and double until accurate
    width = np.sqrt(2 * p.lambda_s * 8.617333262e-5 * p.T) / HARTREE_EV / alpha
    deg = int(min(max_degree, max(40, 8.0 / width)))
    xs = np.linspace(-1, 1, 20001)
    fx = f(xs)
    while True:
        deg += deg % 2
        c = C.chebinterpolate(f, deg)
        err = float(np.max(np.abs(C.chebval(xs, c) - fx)))
        if err < tol or deg >= max_degree:
            return f, c[0::2], deg, rho_max, err
        deg = min(max_degree, 2 * deg)


def _qsvt_prob(H_sec, chi_sec, e_c, alpha, angles, seed=0):
    """||P(A) chi||^2/||chi||^2 from a PennyLane QSVT circuit, A = (H - e_c)/alpha in a random basis."""
    dim = H_sec.shape[0]
    rng = np.random.default_rng(seed)
    Q, _ = np.linalg.qr(rng.standard_normal((dim, dim)) + 1j * rng.standard_normal((dim, dim)))
    Hb = Q.conj().T @ H_sec @ Q
    chi = Q.conj().T @ chi_sec
    chi = chi / np.linalg.norm(chi)
    m = max(1, int(np.ceil(np.log2(dim))))
    pad = 2 ** m
    A = np.eye(pad, dtype=complex) * 0.999
    A[:dim, :dim] = (Hb - e_c * np.eye(dim)) / alpha
    ea, va = np.linalg.eigh(A)
    Sq = (va * np.sqrt(np.clip(1 - ea ** 2, 0, None))) @ va.conj().T
    UA = np.block([[A, Sq], [Sq, -A]])
    assert np.allclose(UA @ UA.conj().T, np.eye(2 * pad), atol=1e-9)
    phi = np.zeros(pad, dtype=complex)
    phi[:dim] = chi
    sys_w = list(range(2, 2 + m))
    anc, ctl = 1, 0
    dev = qml.device("default.qubit", wires=2 + m)

    @qml.qnode(dev)
    def circuit():
        qml.StatePrep(phi, wires=sys_w)
        qml.Hadamard(ctl)
        for sign, cv in ((1, 0), (-1, 1)):
            proj = [qml.PCPhase(sign * a, dim=pad, wires=[anc] + sys_w) for a in angles]
            qml.ctrl(qml.QSVT(qml.QubitUnitary(UA, wires=[anc] + sys_w), proj), control=ctl, control_values=cv)
        qml.Hadamard(ctl)
        return qml.probs(wires=[ctl, anc])

    return float(circuit()[0]), 2 + m


def filtered_isc(asp, states, p: MLJParams = MLJParams(), run_quantum=True, tol=2e-3, max_degree=8000):
    """Gap-resolved k_ISC in the active space: QSVT estimate, exact golden rule, and per-channel data."""
    nq = 2 * asp["ncas"]
    ne = asp["nelec"]
    Hs = states["Hs"]
    Hsoc = qubit_onebody(soc=asp["soc"], n=asp["ncas"])
    Hsoc = qml.simplify(Hsoc)
    coeffs, ops = Hsoc.terms()
    Hsoc_m = qml.dot(np.real(coeffs), ops).sparse_matrix(wire_order=range(nq))
    chi = Hsoc_m @ states["S1"]
    E0 = states["E0"]
    E_S1 = E0 + states["E_S"][1] / HARTREE_EV

    sectors = {}
    lo, hi = np.inf, -np.inf
    for sz2 in (0, 2, -2):
        idx = _sector(nq, ne, sz2)
        Hd = Hs[idx][:, idx].toarray()
        ev = np.linalg.eigvalsh(Hd)
        lo, hi = min(lo, ev[0]), max(hi, ev[-1])
        sectors[sz2] = (idx, Hd)
    e_c = lo - 1e-3                      # all spectra at x >= 0: one even polynomial suffices
    alpha = 1.01 * (hi - e_c)
    f, ceven, deg, rho_max, approx_err = filter_poly(alpha, e_c, E_S1, p, tol=tol, max_degree=max_degree)
    angles = None
    if run_quantum:
        phases, fpi_err, _ = qsp_phases_fpi(ceven, deg)
        angles = to_pennylane_qsvt(phases)

    pref = 2 * np.pi / HBAR_EVS * HARTREE_EV ** 2      # chi in Ha -> V^2 in eV^2
    out = {"alpha_Ha": float(alpha), "degree": int(deg), "poly_err": approx_err, "rho_max_per_eV": rho_max,
           "E_S1_eV": float(states["E_S"][1]), "sectors": {}, "channels": []}
    k_q = k_exact = k_poly = 0.0
    for sz2, (idx, Hd) in sectors.items():
        chi_s = chi[idx]
        n2 = float(np.vdot(chi_s, chi_s).real)
        e, v = np.linalg.eigh(Hd)
        ov = np.abs(v.conj().T @ chi_s) ** 2
        dE = (e - E_S1) * HARTREE_EV
        rho = _rho_vec(dE, p)
        k_exact += pref * float(np.sum(ov * rho))
        k_poly += pref * n2 * rho_max * float(np.sum(ov / max(n2, 1e-300) *
                                                      np.polynomial.chebyshev.chebval((e - e_c) / alpha,
                                                                                      _full(ceven)) ** 2)) / 0.98 ** 2
        sec = {"chi_norm2_Ha2": n2, "dim": int(len(idx))}
        if run_quantum and n2 > 1e-20:
            prob, nw = _qsvt_prob(Hd, chi_s, e_c, alpha, angles)
            sec.update(p_success=prob, wires=nw)
            k_q += pref * n2 * rho_max * prob / 0.98 ** 2
        out["sectors"][int(sz2)] = sec
        # channel list (sector 0 shown once: the M=0 triplets; weights are per M)
        if sz2 == 0:
            for j in np.argsort(dE):
                if ov[j] > 1e-16 * max(n2, 1e-300) and n2 > 0:
                    out["channels"].append({"dE_eV": float(dE[j]), "V2_M0_cm2": float(ov[j]) * CM1_PER_HARTREE ** 2,
                                            "k_M0": pref * float(ov[j]) * float(rho[j])})
    out.update(k_isc_qsvt=k_q if run_quantum else None, k_isc_exact=k_exact, k_isc_poly=k_poly)
    # the gap-free proxy of arXiv:2512.15889 for comparison: |<T1|H_SO|S1>|^2 summed over M
    T1 = states["T"][0]
    V2_T1 = sum(abs(np.vdot(T1[M], chi)) ** 2 for M in (-1, 0, 1))
    out["proxy_T1_V_cm"] = float(np.sqrt(V2_T1) * CM1_PER_HARTREE)
    out["dE_T1_eV"] = float(states["E_T"][0] - states["E_S"][1])
    return out


def _full(ceven):
    c = np.zeros(2 * len(ceven) - 1)
    c[0::2] = ceven
    return c


def resource_estimate(lambda_Ha, toffoli_per_walk, p: MLJParams = MLJParams(), e_th_Ha=0.0):
    """Toffolis per shot of the filtered observable: polynomial degree from the empirical fit of
    arXiv:2512.15889 (d = 4.7571 lambda'/Delta + 321.2, lambda' = lambda + |E_th|) with Delta the
    filter width, times the cost of one walk step."""
    delta = np.sqrt(2 * p.lambda_s * 8.617333262e-5 * p.T) / HARTREE_EV
    lam = lambda_Ha + abs(e_th_Ha)
    d = 4.7571 * lam / delta + 321.2051
    return {"filter_width_eV": float(delta * HARTREE_EV), "degree": float(d), "toffoli_per_shot": float(d * toffoli_per_walk)}
