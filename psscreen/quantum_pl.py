"""Backend cuántico en PennyLane siguiendo arXiv:2512.15889 (Xanadu), ejecutado en simulador.

Dos algoritmos del paper, sobre un espacio activo pequeño (CAS(n,n) -> 2n qubits, Jordan-Wigner):

1) Proyección por umbral para absorción acumulada en una ventana de energías.
   Paper: |phi> = mu|Psi0>/||mu Psi0||, block-encoding (THC + qubitización) y QSP generalizado
   con dos escalones de Heaviside; P(éxito) = peso dipolar del espectro dentro de la ventana.
   Aquí: block-encoding por inmersión (qml.BlockEncode) del Hamiltoniano restringido al sector
   N fijo, S=0 (base de funciones de espín adaptadas, equivalente a CSFs), y QSVT (qml.QSVT +
   qml.PCPhase) con UN polinomio par que aproxima la ventana completa sobre x=(H-E_c)/alpha. La
   parte real del bloque se obtiene con una LCU de dos QSVT (ángulos +phi y -phi). Ángulos de grado
   alto por iteración de punto fijo (qsp_angles.py). Mismo observable que el paper.

2) Proxy de ISC por evolución temporal.
   Paper: k(t) = |<T1,M| exp(-i H_SOC t) |S1>|^2  ~ t^2 |<T1,M|H_SOC|S1>|^2 a tiempos cortos.
   Aquí: H_SOC monoelectrónico en espín-orbitales (integrales calibradas de soc.py) -> qubits,
   evolución con qml.TrotterProduct, estados S1 y T1(M=-1,0,+1) del Hamiltoniano activo exacto
   cargados con qml.StatePrep (el paper los prepara con sum-of-Slaters).

Además (estrategia de Atchade-Adelomou & Coronas, Electronics 2025): tapering Z2 del Hamiltoniano
de qubits y estimación de recursos FTQC con qml.resource para espacios activos mayores.
"""
from __future__ import annotations
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import pennylane as qml
from pyscf import mcscf, ao2mo
from .qsp_angles import qsp_phases_fpi, to_pennylane_qsvt
from .soc import soc_mo_integrals, CM1_PER_HARTREE

HARTREE_EV = 27.211386245988
EV_NM = 1239.84198


# ----------------------------------------------------------------------------- espacio activo
def active_space(mf, ncas, nelecas, basis):
    mc = mcscf.CASCI(mf, ncas, nelecas)
    h1, ecore = mc.get_h1eff()
    h2 = ao2mo.restore(1, mc.get_h2eff(), ncas)
    ncore = mc.ncore
    idx = list(range(ncore, ncore + ncas))
    C = mf.mo_coeff[:, idx]
    mol = mf.mol
    with mol.with_common_orig((0, 0, 0)):
        r = mol.intor_symmetric("int1e_r", comp=3)
    dip = -np.einsum("xpq,pi,qj->xij", r, C, C)   # operador dipolar electrónico (carga -1)
    h_soc_mo, scales = soc_mo_integrals(mol, mf.mo_coeff, basis)
    h_soc = h_soc_mo[:, idx][:, :, idx]
    return dict(h1=h1, h2=h2, ecore=float(ecore), dip=dip, soc=h_soc, ncas=ncas,
                nelec=int(np.sum(nelecas)) if not np.isscalar(nelecas) else int(nelecas),
                orb_energies_eV=(mf.mo_energy[idx] * HARTREE_EV).tolist())


def qubit_hamiltonian(asp):
    op = qml.qchem.fermionic_observable(np.array([asp["ecore"]]), asp["h1"], asp["h2"].transpose(0, 2, 3, 1))
    return qml.jordan_wigner(op)


def qubit_onebody(mat_spin_summed=None, soc=None, n=None):
    """Operador de un cuerpo en espín-orbitales (wire 2p = p alfa, 2p+1 = p beta)."""
    fs = qml.fermi.FermiSentence({})
    if mat_spin_summed is not None:
        n = mat_spin_summed.shape[0]
        for p in range(n):
            for q in range(n):
                if abs(mat_spin_summed[p, q]) > 1e-12:
                    for s in (0, 1):
                        fs += mat_spin_summed[p, q] * qml.FermiWord({(0, 2 * p + s): "+", (1, 2 * q + s): "-"})
    if soc is not None:
        n = soc.shape[1]
        pauli = [np.array([[0, 1], [1, 0]]) / 2, np.array([[0, -1j], [1j, 0]]) / 2, np.array([[1, 0], [0, -1]]) / 2]
        hs = sum(np.einsum("pq,st->stpq", soc[k], pauli[k]) for k in range(3))   # h^{st}_pq
        for s in (0, 1):
            for t in (0, 1):
                for p in range(n):
                    for q in range(n):
                        c = hs[s, t, p, q]
                        if abs(c) > 1e-14:
                            fs += c * qml.FermiWord({(0, 2 * p + s): "+", (1, 2 * q + t): "-"})
    return qml.jordan_wigner(fs)


def _sector(nq, nelec, sz2):
    """Índices de la base computacional con N electrones y 2*Sz = sz2 (wire 0 = bit más significativo)."""
    out = []
    for i in range(2 ** nq):
        bits = [(i >> (nq - 1 - w)) & 1 for w in range(nq)]
        na, nb = sum(bits[0::2]), sum(bits[1::2])
        if na + nb == nelec and na - nb == sz2:
            out.append(i)
    return np.array(out)


# ----------------------------------------------------------------------------- estados exactos
def eigenstates(asp, H=None):
    nq = 2 * asp["ncas"]; ne = asp["nelec"]
    H = H if H is not None else qubit_hamiltonian(asp)
    Hs = H.sparse_matrix(wire_order=range(nq)).tocsr()
    S2 = qml.qchem.spin2(ne, nq).sparse_matrix(wire_order=range(nq)).tocsr()
    res = {}
    for sz2 in (0, 2, -2):
        idx = _sector(nq, ne, sz2)
        Hd = Hs[idx][:, idx].toarray(); Sd = S2[idx][:, idx].toarray()
        e, v = np.linalg.eigh(Hd)
        s2 = np.einsum("ik,ij,jk->k", v.conj(), Sd, v).real
        res[sz2] = (idx, e, v, s2)
    idx0, e0, v0, s20 = res[0]
    sing = [k for k in range(len(e0)) if abs(s20[k]) < 1e-3]
    trip0 = [k for k in range(len(e0)) if abs(s20[k] - 2) < 1e-3]

    def full(idx, vec):
        f = np.zeros(2 ** nq, dtype=complex); f[idx] = vec; return f
    E0 = e0[sing[0]]
    states = {"E0": E0, "Psi0": full(idx0, v0[:, sing[0]]), "S1": full(idx0, v0[:, sing[1]]),
              "E_S": (e0[sing] - E0) * HARTREE_EV, "E_T": (e0[trip0] - E0) * HARTREE_EV,
              "T": []}
    ntrip = min(3, len(trip0))
    for n in range(ntrip):
        states["T"].append({0: full(idx0, v0[:, trip0[n]])})
    for sz2, M in ((2, 1), (-2, -1)):
        idx, e, v, s2 = res[sz2]
        tr = [j for j in range(len(e)) if abs(s2[j] - 2) < 1e-3]
        for n in range(ntrip):
            assert abs(e[tr[n]] - e0[trip0[n]]) < 1e-6, "T_n no degenerado entre M (casi-degeneración?)"
            states["T"][n][M] = full(idx, v[:, tr[n]])
    states["T1"] = states["T"][0]
    states["Hs"], states["S2"], states["sector0"] = Hs, S2, (idx0, e0, v0, s20)
    return states


# ----------------------------------------------------------------------------- 1) absorción QSVT
def window_poly(x_lo, x_hi, sigma, degree, amp=0.98):
    """Polinomio par en Chebyshev ~ amp * ventana suave en |x| in [x_lo, x_hi] (x_lo puede ser 0)."""
    from numpy.polynomial import chebyshev as C
    from scipy.special import erf

    def f(x):
        a = np.abs(x)
        return amp * 0.25 * (1 + erf((a - x_lo) / sigma)) * (1 + erf((x_hi - a) / sigma)) if x_lo > 0 else \
            amp * 0.5 * (1 + erf((x_hi - a) / sigma))
    c = C.chebinterpolate(f, degree)
    return f, c[0::2]


def threshold_projection(asp, states, e_lo_ev, e_hi_ev, edge_ev=0.05, max_degree=2000, shots=None):
    """Peso dipolar en la ventana [e_lo, e_hi] (energía de excitación) con QSVT en PennyLane."""
    idx0, e0, v0, s20 = states["sector0"]
    sing = [k for k in range(len(e0)) if abs(s20[k]) < 1e-3]
    B = v0[:, sing]                                          # base S=0 (autoestados = CSF-like)
    Hsing = np.diag(e0[sing])                               # H en esa base (classically diag: ver nota)
    # para no "hacer trampa" con la base propia, rotamos a una base S=0 genérica (aleatoria fija)
    rng = np.random.default_rng(0)
    Q, _ = np.linalg.qr(rng.standard_normal((len(sing), len(sing))))
    Hb = Q.T @ Hsing @ Q
    nq = 2 * asp["ncas"]
    mus = []
    for a in range(3):
        D = qubit_onebody(mat_spin_summed=asp["dip"][a]).sparse_matrix(wire_order=range(nq)).tocsr()
        Dd = D[idx0][:, idx0].toarray()
        mus.append(Q.T @ (B.conj().T @ Dd @ B).real @ Q)
    psi0 = Q.T @ np.eye(len(sing))[:, 0]
    E0 = e0[sing[0]]
    e_c = E0 + 0.5 * (e_lo_ev + e_hi_ev) / HARTREE_EV
    evals = np.linalg.eigvalsh(Hb)
    alpha = 1.01 * np.max(np.abs(evals - e_c))
    w_half = 0.5 * (e_hi_ev - e_lo_ev) / HARTREE_EV / alpha
    sigma = edge_ev / HARTREE_EV / alpha
    if 6.0 / sigma > max_degree:        # el grado disponible no resuelve ese borde: se ensancha
        sigma = 6.0 / max_degree
        edge_ev = sigma * alpha * HARTREE_EV
    degree = int(max(20, 6.0 / sigma))
    degree += degree % 2
    f, ceven = window_poly(0.0, w_half, sigma, degree)
    phases, err, _ = qsp_phases_fpi(ceven, degree)
    ang = to_pennylane_qsvt(phases)
    dim = len(sing); m = int(np.ceil(np.log2(dim))) if dim > 1 else 1; pad = 2 ** m
    A = np.eye(pad) * 0.999                                   # relleno fuera de la ventana
    A[:dim, :dim] = (Hb - e_c * np.eye(dim)) / alpha
    sys_w = list(range(2, 2 + m)); anc, ctl = 1, 0
    # Unitario de block-encoding U = [[A, sqrt(1-A^2)], [sqrt(1-A^2), -A]] (el mismo que define
    # qml.BlockEncode). Se construye con la descomposición espectral de A porque scipy.sqrtm, que usa
    # qml.BlockEncode, pierde precisión en matrices > ~100x100 (error de bloque ~0.3 medido aquí).
    ea, va = np.linalg.eigh(A)
    Sq = (va * np.sqrt(np.clip(1 - ea ** 2, 0, None))) @ va.T
    UA = np.block([[A, Sq], [Sq, -A]])
    assert np.allclose(UA @ UA.T, np.eye(2 * pad), atol=1e-10)
    dev = qml.device("default.qubit", wires=2 + m, shots=shots)
    out = {"alpha_Ha": alpha, "degree": int(degree), "fpi_err": float(err), "n_qubits_sim": 2 + m,
           "n_qubits_fock": nq, "dim_S0": dim, "window_eV": (e_lo_ev, e_hi_ev), "edge_eV": edge_ev}

    @qml.qnode(dev)
    def circuit(phi):
        qml.StatePrep(phi, wires=sys_w)
        qml.Hadamard(ctl)
        for sign, cv in ((1, 0), (-1, 1)):
            proj = [qml.PCPhase(sign * a, dim=pad, wires=[anc] + sys_w) for a in ang]
            qml.ctrl(qml.QSVT(qml.QubitUnitary(UA, wires=[anc] + sys_w), proj), control=ctl, control_values=cv)
        qml.Hadamard(ctl)
        return qml.probs(wires=[ctl, anc])

    xs = (evals - e_c) / alpha
    ev_vecs = np.linalg.eigh(Hb)[1]
    tot_w, q_w, ex_smooth, ex_sharp = 0.0, 0.0, 0.0, 0.0
    for a in range(3):
        phi = mus[a] @ psi0 - (psi0 @ mus[a] @ psi0) * psi0        # parte de transición
        w = float(phi @ phi)
        if w < 1e-14:
            continue
        phi_n = np.zeros(pad); phi_n[:dim] = phi / np.sqrt(w)
        p = circuit(phi_n)[0]                                   # P(ctl=0, anc=0) = ||f(A) phi||^2
        ov = (ev_vecs.T @ phi_n[:dim]) ** 2
        tot_w += w; q_w += w * p
        ex_smooth += w * float(np.sum(ov * f(xs) ** 2))
        Eexc = (evals - E0) * HARTREE_EV
        ex_sharp += w * float(np.sum(ov[(Eexc >= e_lo_ev) & (Eexc <= e_hi_ev)]) * 0.98 ** 2)
    out["edge_eV"] = float(edge_ev)
    out.update(weight_quantum=q_w / tot_w / 0.98 ** 2, weight_exact_smooth=ex_smooth / tot_w / 0.98 ** 2,
               weight_exact_sharp=ex_sharp / tot_w / 0.98 ** 2, dipole_norm2=tot_w)
    return out


# ----------------------------------------------------------------------------- 2) proxy ISC
def isc_proxy(asp, states, times=(20.0, 40.0), trotter_steps=4, triplet=0):
    nq = 2 * asp["ncas"]
    Hsoc = qubit_onebody(soc=asp["soc"], n=asp["ncas"])
    Hsoc = qml.simplify(Hsoc)
    coeffs, ops = Hsoc.terms()
    assert np.allclose(np.imag(coeffs), 0, atol=1e-12), "H_SOC en qubits debe ser hermítico (coef. reales)"
    Hsoc = qml.dot(np.real(coeffs), ops)
    dev = qml.device("default.qubit", wires=nq)

    @qml.qnode(dev)
    def evolve(t):
        qml.StatePrep(states["S1"], wires=range(nq))
        qml.TrotterProduct(Hsoc, time=t, n=trotter_steps, order=2)
        return qml.state()

    # referencia clásica exacta
    Hm = Hsoc.sparse_matrix(wire_order=range(nq))
    T = states["T"][triplet]
    V2 = sum(abs(np.vdot(T[M], Hm @ states["S1"])) ** 2 for M in (-1, 0, 1))
    res = {"triplet": triplet + 1, "V_exact_cm": float(np.sqrt(V2) * CM1_PER_HARTREE), "n_pauli_soc": len(ops), "t": [], "k": [], "V_proxy_cm": []}
    for t in times:
        psi = evolve(t)
        k = sum(abs(np.vdot(T[M], psi)) ** 2 for M in (-1, 0, 1))
        res["t"].append(t); res["k"].append(float(k)); res["V_proxy_cm"].append(float(np.sqrt(k) / t * CM1_PER_HARTREE))
    return res


# ----------------------------------------------------------------------------- tapering y recursos
def taper_report(asp, H=None):
    nq = 2 * asp["ncas"]
    H = H if H is not None else qubit_hamiltonian(asp)
    gens = qml.qchem.symmetry_generators(H)
    paulixs = qml.qchem.paulix_ops(gens, nq)
    sector = qml.qchem.optimal_sector(H, gens, asp["nelec"])
    Ht = qml.taper(H, gens, paulixs, sector)
    return {"qubits": nq, "qubits_tapered": len(Ht.wires), "n_generators": len(gens),
            "pauli_terms": len(H.terms()[0]), "pauli_terms_tapered": len(Ht.terms()[0])}


def ftqc_resources(asp, error=0.0016):
    """QPE por qubitización con doble factorización (qml.resource) para este espacio activo."""
    from pennylane.estimator.qpe_resources.second_quantization import DoubleFactorization
    # PennyLane espera las integrales en el mismo orden que electron_integrals (como fermionic_observable)
    df = DoubleFactorization(asp["h1"], asp["h2"].transpose(0, 2, 3, 1), error=error)
    return {"logical_qubits": int(df.qubits), "toffoli_gates": float(df.gates), "lambda": float(df.lamb)}
