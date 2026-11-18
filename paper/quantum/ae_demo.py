"""Amplitude estimation on the gap-resolved ISC circuit (PennyLane, simulator).

Builds the unitary A of the filter circuit (state preparation of chi, LCU of two QSVT sequences) for one M_S
sector of the BODIPY-core CAS(4,4) problem and the Grover operator Q = -A S_0 A^dagger S_good. Then compares, at
equal numbers of calls to A (or A^dagger):

  1. direct sampling of the success flag (exact binomial statistics);
  2. canonical QAE (QPE on Q, median of 16 runs). With --qpe the QPE circuit is simulated with
     qml.QuantumPhaseEstimation (t = 4..13); the analytic QPE outcome distribution is used for the ensemble below.
     The true p puts 2^t theta / pi = 3.07 at t = 12, close to a grid point, so the circuit result is a
     favourable case; we therefore also report canonical QAE averaged over an ensemble of nearby amplitudes
     theta' = theta (1 + u), u ~ U(-0.3, 0.3);
  3. maximum-likelihood amplitude estimation (MLAE, Suzuki et al. 2020) with the exponential schedule
     m = 0, 1, 2, 4, ..., 2^K and N_SHOTS shots per depth; no QPE, no controlled Q.

The good-state probability after m Grover iterations is checked against sin^2((2m+1) theta) on the actual Q, which
justifies simulating QAE and MLAE from theta. Usage (repository root):
    python paper/quantum/ae_demo.py [--qpe] [NAME NCAS SZ2]
Writes paper/quantum/results/ae_demo.json and paper/figures/ae_demo.pdf.
"""
import json
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
import numpy as np  # noqa: E402
from scipy.stats import binom  # noqa: E402

N_REP = 16          # median of 16 canonical-QAE runs (confidence 0.95, Brassard et al. Thm 12 + Hoeffding)
N_SHOTS = 100       # MLAE shots per Grover depth
N_TRIALS = 2000     # Monte Carlo trials for MLAE and for the median-of-16 canonical QAE
CONF = 0.95


def build(name="BODIPY_core", ncas=4, sz2_sel=2):
    import pennylane as qml
    from pyscf import dft, gto

    from psscreen.backends import _ecp
    from psscreen.geometry import read_xyz, to_pyscf_atom
    from psscreen.isc_filtered import filter_poly
    from psscreen.qsp_angles import qsp_phases_fpi, to_pennylane_qsvt
    from psscreen.quantum_pl import (HARTREE_EV, _sector, active_space, eigenstates, qubit_hamiltonian,
                                     qubit_onebody)
    from psscreen.rates import MLJParams

    sym, xyz = read_xyz(os.path.join(ROOT, "validation", "runs", f"{name}.xyz"))
    atom = to_pyscf_atom(sym, xyz)
    mol = gto.M(atom=atom, basis="def2-svp", ecp=_ecp(atom, "def2-svp"), verbose=0)
    mf = dft.RKS(mol, xc="b3lyp").density_fit()
    mf.grids.level = 2
    mf.kernel()
    asp = active_space(mf, ncas, ncas, "def2-svp")
    st = eigenstates(asp, qubit_hamiltonian(asp))
    nq, ne = 2 * ncas, asp["nelec"]
    Hsoc = qml.simplify(qubit_onebody(soc=asp["soc"], n=ncas))
    co, ops = Hsoc.terms()
    chi = qml.dot(np.real(co), ops).sparse_matrix(wire_order=range(nq)) @ st["S1"]
    lo, hi, blocks = np.inf, -np.inf, {}
    for sz2 in (0, 2, -2):
        idx = _sector(nq, ne, sz2)
        Hd = st["Hs"][idx][:, idx].toarray()
        ev = np.linalg.eigvalsh(Hd)
        lo, hi = min(lo, ev[0]), max(hi, ev[-1])
        blocks[sz2] = (idx, Hd)
    e_c = lo - 1e-3
    alpha = 1.01 * (hi - e_c)
    E_S1 = st["E0"] + st["E_S"][1] / HARTREE_EV
    f, ceven, deg, rho_max, err = filter_poly(alpha, e_c, E_S1, MLJParams())
    angles = to_pennylane_qsvt(qsp_phases_fpi(ceven, deg)[0])

    idx, Hd = blocks[sz2_sel]
    dim = Hd.shape[0]
    rng = np.random.default_rng(0)
    Qr, _ = np.linalg.qr(rng.standard_normal((dim, dim)) + 1j * rng.standard_normal((dim, dim)))
    Hb = Qr.conj().T @ Hd @ Qr
    c = Qr.conj().T @ chi[idx]
    c = c / np.linalg.norm(c)
    m = max(1, int(np.ceil(np.log2(dim))))
    pad = 2 ** m
    A_ = np.eye(pad, dtype=complex) * 0.999
    A_[:dim, :dim] = (Hb - e_c * np.eye(dim)) / alpha
    ea, va = np.linalg.eigh(A_)
    Sq = (va * np.sqrt(np.clip(1 - ea ** 2, 0, None))) @ va.conj().T
    UA = np.block([[A_, Sq], [Sq, -A_]])
    phi = np.zeros(pad, dtype=complex)
    phi[:dim] = c
    sys_w = list(range(2, 2 + m))
    wires = list(range(2 + m))

    def A_circuit():
        qml.StatePrep(phi, wires=sys_w)
        qml.Hadamard(0)
        for sign, cv in ((1, 0), (-1, 1)):
            proj = [qml.PCPhase(sign * a, dim=pad, wires=[1] + sys_w) for a in angles]
            qml.ctrl(qml.QSVT(qml.QubitUnitary(UA, wires=[1] + sys_w), proj), control=0, control_values=cv)
        qml.Hadamard(0)

    A = qml.matrix(A_circuit, wire_order=wires)()
    n = len(wires)
    good = np.array([((i >> (n - 1)) & 1) == 0 and ((i >> (n - 2)) & 1) == 0 for i in range(2 ** n)])
    p_true = float(np.sum(np.abs(A[:, 0][good]) ** 2))
    S_good = np.diag(np.where(good, -1.0, 1.0))
    S_0 = np.eye(2 ** n)
    S_0[0, 0] = -1.0
    Q = -A @ S_0 @ A.conj().T @ S_good
    return dict(A=A, Q=Q, good=good, n=n, wires=wires, p=p_true, degree=deg, dim=dim, name=name, ncas=ncas,
                sz2=sz2_sel)


def check_grover_law(pr, ms=(0, 1, 2, 5, 10, 50, 100, 200, 331, 500)):
    theta = math.asin(math.sqrt(pr["p"]))
    psi = pr["A"][:, 0].copy()
    out, m_prev = [], 0
    for m in ms:
        for _ in range(m - m_prev):
            psi = pr["Q"] @ psi
        m_prev = m
        pg = float(np.sum(np.abs(psi[pr["good"]]) ** 2))
        out.append({"m": m, "p_good": pg, "sin2": math.sin((2 * m + 1) * theta) ** 2})
    return out


def qpe_dist(theta, t):
    """Outcome distribution of canonical QAE with t evaluation qubits (Brassard et al. 2002)."""
    T = 2 ** t
    y = np.arange(T)
    pr = np.zeros(T)
    for s in (1, -1):
        d = y / T - s * theta / math.pi
        num = np.sin(T * math.pi * d) ** 2
        den = (T ** 2) * np.sin(math.pi * d) ** 2
        with np.errstate(invalid="ignore", divide="ignore"):
            v = np.where(den < 1e-300, 1.0, num / den)
        pr += 0.5 * v
    return pr / pr.sum()


def canonical_qae_err(theta, t, rng, n_trials=N_TRIALS):
    p = math.sin(theta) ** 2
    probs = qpe_dist(theta, t)
    p_hat = np.sin(np.pi * np.arange(2 ** t) / 2 ** t) ** 2
    draws = rng.choice(p_hat, size=(n_trials, N_REP), p=probs)
    return np.abs(np.median(draws, axis=1) - p) / p


def mlae_err(theta, K, rng, n_trials=N_TRIALS, shots=N_SHOTS):
    """MLAE with depths m = 0, 1, 2, 4, ..., 2^K; hierarchical grid maximum likelihood."""
    ms = [0] + [2 ** k for k in range(K + 1)]
    p = math.sin(theta) ** 2
    hits = [rng.binomial(shots, math.sin((2 * m + 1) * theta) ** 2, size=n_trials) for m in ms]
    est = np.empty(n_trials)
    for i in range(n_trials):
        lo, hi = 1e-9, math.pi / 2 - 1e-9
        fisher = 0.0
        for j, m in enumerate(ms):
            fisher += shots * (2 * (2 * m + 1)) ** 2
            g = np.linspace(lo, hi, 12001)
            ll = np.zeros_like(g)
            for jj in range(j + 1):
                q = np.clip(np.sin((2 * ms[jj] + 1) * g) ** 2, 1e-15, 1 - 1e-15)
                ll += hits[jj][i] * np.log(q) + (shots - hits[jj][i]) * np.log(1 - q)
            th = g[np.argmax(ll)]
            w = 30 / math.sqrt(fisher)
            lo, hi = max(1e-9, th - w), min(math.pi / 2 - 1e-9, th + w)
        est[i] = th
    calls = sum(shots * (2 * m + 1) for m in ms)
    return np.abs(np.sin(est) ** 2 - p) / p, calls


def sampling_err_quantile(p, calls, q):
    k = np.arange(0, min(calls, int(50 * calls * p + 200)) + 1)
    pk = binom.pmf(k, calls, p)
    errs = np.abs(k / calls - p) / p
    order = np.argsort(errs)
    cdf = np.cumsum(pk[order])
    return float(errs[order][min(np.searchsorted(cdf, q), len(order) - 1)])


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    run_qpe = "--qpe" in sys.argv
    name = args[0] if args else "BODIPY_core"
    ncas = int(args[1]) if len(args) > 1 else 4
    sz2 = int(args[2]) if len(args) > 2 else 2
    out_path = os.path.join(ROOT, "paper", "quantum", "results", "ae_demo.json")
    old = json.load(open(out_path)) if os.path.exists(out_path) else {}

    pr = build(name, ncas, sz2)
    p = pr["p"]
    theta = math.asin(math.sqrt(p))
    print(f"{name} CAS({ncas},{ncas}) 2Sz={sz2}: wires {pr['n']}, degree {pr['degree']}, p = {p:.4e}", flush=True)
    law = check_grover_law(pr)
    dev_max = max(abs(r["p_good"] - r["sin2"]) for r in law)
    print(f"Grover law on the actual Q: max |p_good - sin^2((2m+1)theta)| = {dev_max:.1e}", flush=True)

    circuit_qpe = old.get("circuit_qpe", old.get("results"))
    if run_qpe:
        circuit_qpe = run_circuit_qpe(pr)

    rng = np.random.default_rng(7)
    ens = theta * (1 + rng.uniform(-0.3, 0.3, 64))
    canon = []
    for t in range(6, 17):
        M = 2 ** t - 1
        calls = N_REP * (2 * M + 1)
        e_true = canonical_qae_err(theta, t, rng)
        e_ens = np.concatenate([canonical_qae_err(th, t, rng, 200) for th in ens])
        canon.append({"t": t, "calls": calls, "q95_true": float(np.quantile(e_true, CONF)),
                      "q95_ensemble": float(np.quantile(e_ens, CONF)),
                      "median_ensemble": float(np.median(e_ens)),
                      "q95_sampling": sampling_err_quantile(p, calls, CONF)})
        print(f"canonical t={t:2d} calls={calls:8d} q95 rel err: true-theta {canon[-1]['q95_true']:.3f}, "
              f"ensemble {canon[-1]['q95_ensemble']:.3f}; sampling {canon[-1]['q95_sampling']:.3f}", flush=True)

    mlae = []
    for K in range(4, 14):
        e_true, calls = mlae_err(theta, K, rng, 400)
        e_ens = np.concatenate([mlae_err(th, K, rng, 25)[0] for th in ens[:16]])
        rec = {"K": K, "calls": calls, "q95_true": float(np.quantile(e_true, CONF)),
               "q95_ensemble": float(np.quantile(e_ens, CONF)), "median_ensemble": float(np.median(e_ens)),
               "q95_sampling": sampling_err_quantile(p, calls, CONF)}
        rec["c_eff"] = rec["q95_ensemble"] * math.sqrt(p) * calls       # eps_rel ~ c / (sqrt(p) N_calls)
        mlae.append(rec)
        print(f"MLAE K={K:2d} calls={calls:8d} q95 rel err: true {rec['q95_true']:.3f}, ensemble "
              f"{rec['q95_ensemble']:.3f} (c = {rec['c_eff']:.1f}); sampling {rec['q95_sampling']:.3f}", flush=True)

    # empirical constants in the Heisenberg regime (95th-percentile relative error below 0.5)
    c_mlae = float(np.median([r["c_eff"] for r in mlae if r["q95_ensemble"] < 0.5]))
    c_canon = float(np.median([r["q95_ensemble"] * math.sqrt(p) * r["calls"] for r in canon
                               if r["q95_ensemble"] < 0.5]))
    out = {"name": name, "ncas": ncas, "sz2": sz2, "degree": pr["degree"], "p": p, "confidence": CONF,
           "grover_law_check": law, "grover_law_max_dev": dev_max, "circuit_qpe": circuit_qpe,
           "canonical_analytic": canon, "mlae": mlae, "mlae_shots_per_depth": N_SHOTS,
           "c_mlae_q95": c_mlae, "c_canonical_q95": c_canon,
           "heisenberg_threshold_iterations": math.pi / math.sqrt(p)}
    print(f"empirical constants (eps_rel * sqrt(p) * calls at the 95th percentile): MLAE {c_mlae:.1f}, "
          f"canonical median-of-16 {c_canon:.1f}", flush=True)
    json.dump(out, open(out_path, "w"), indent=1)
    fig(out)


def run_circuit_qpe(pr):
    import pennylane as qml
    p = pr["p"]
    res = []
    for t in range(4, 14):
        est = list(range(pr["n"], pr["n"] + t))
        dev = qml.device("default.qubit", wires=pr["n"] + t)

        @qml.qnode(dev)
        def qae():
            qml.QubitUnitary(pr["A"], wires=pr["wires"])
            qml.QuantumPhaseEstimation(qml.QubitUnitary(pr["Q"], wires=pr["wires"]), estimation_wires=est)
            return qml.probs(wires=est)

        probs = np.clip(np.asarray(qae()), 0, None)
        probs /= probs.sum()
        p_hat = np.sin(np.pi * np.arange(2 ** t) / 2 ** t) ** 2
        calls = N_REP * (2 * (2 ** t - 1) + 1)
        draws = np.random.default_rng(1).choice(p_hat, size=(4000, N_REP), p=probs)
        e = np.abs(np.median(draws, axis=1) - p)
        res.append({"t": t, "calls": calls, "median_abs_err_ae": float(np.median(e)),
                    "q90_abs_err_ae": float(np.quantile(e, 0.9))})
    return res


def fig(out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    p = out["p"]
    fig_, ax = plt.subplots(figsize=(3.4, 3.6))
    c = out["canonical_analytic"]
    m = out["mlae"]
    ax.loglog([r["calls"] for r in c], [max(r["q95_sampling"], 1e-4) for r in c], "s-", color="#c05621", ms=3,
              label="direct sampling")
    ax.loglog([r["calls"] for r in c], [max(r["q95_ensemble"], 1e-4) for r in c], "o-", color="#2b6cb0", ms=3,
              label="canonical QAE (median of 16)")
    ax.loglog([r["calls"] for r in m], [max(r["q95_ensemble"], 1e-4) for r in m], "^-", color="#2f855a", ms=3,
              label="MLAE")
    if out.get("circuit_qpe"):
        q = out["circuit_qpe"]
        ax.loglog([r["calls"] for r in q], [max(r["q90_abs_err_ae"] / p, 1e-4) for r in q], "x", color="#2b6cb0",
                  ms=4, label="canonical QAE, circuit (this $p$)")
    ax.axhline(0.2, color="0.6", ls=":", lw=0.8)
    ys = [r["q95_sampling"] for r in c] + [r["q95_ensemble"] for r in c + m]
    ax.set_ylim(1e-2, 2 * max(ys))
    ax.set_xlabel("calls to the filter circuit")
    ax.set_ylabel(r"95th-percentile relative error in $k_\mathrm{ISC}$")
    ax.set_title(f"{out['name'].replace('_', ' ')}, CAS({out['ncas']},{out['ncas']}), $p$ = {p:.1e}", fontsize=7)
    ax.legend(frameon=False, fontsize=6, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2)
    fig_.tight_layout()
    fig_.savefig(os.path.join(ROOT, "paper", "figures", "ae_demo.pdf"))


if __name__ == "__main__":
    main()
