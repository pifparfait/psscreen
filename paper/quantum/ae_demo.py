"""Amplitude estimation on the gap-resolved ISC circuit (PennyLane, simulator).

Builds the unitary A of the filter circuit (state preparation of chi, LCU of two QSVT sequences) for one
M_S sector of the BODIPY-core CAS(4,4) problem, forms the Grover operator Q = A S_0 A^dagger S_good, runs
canonical quantum amplitude estimation (qml.QuantumPhaseEstimation on Q) and compares its error with direct
sampling of the success flag at equal numbers of calls to A.

Usage (from the repository root): python paper/quantum/ae_demo.py [NAME] [NCAS] [SZ2]
Writes paper/quantum/results/ae_demo.json and paper/figures/ae_demo.pdf.
"""
import json
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
import numpy as np  # noqa: E402
import pennylane as qml  # noqa: E402
from pyscf import dft, gto  # noqa: E402

from psscreen.backends import _ecp  # noqa: E402
from psscreen.geometry import read_xyz, to_pyscf_atom  # noqa: E402
from psscreen.isc_filtered import filter_poly  # noqa: E402
from psscreen.qsp_angles import qsp_phases_fpi, to_pennylane_qsvt  # noqa: E402
from psscreen.quantum_pl import (HARTREE_EV, _sector, active_space, eigenstates, qubit_hamiltonian,  # noqa: E402
                                 qubit_onebody)
from psscreen.rates import MLJParams  # noqa: E402

name = sys.argv[1] if len(sys.argv) > 1 else "BODIPY_core"
ncas = int(sys.argv[2]) if len(sys.argv) > 2 else 4
sz2_sel = int(sys.argv[3]) if len(sys.argv) > 3 else 2

sym, xyz = read_xyz(os.path.join(ROOT, "validation", "runs", f"{name}.xyz"))
atom = to_pyscf_atom(sym, xyz)
mol = gto.M(atom=atom, basis="def2-svp", ecp=_ecp(atom, "def2-svp"), verbose=0)
mf = dft.RKS(mol, xc="b3lyp").density_fit()
mf.grids.level = 2
mf.kernel()
asp = active_space(mf, ncas, ncas, "def2-svp")
st = eigenstates(asp, qubit_hamiltonian(asp))
nq, ne = 2 * ncas, asp["nelec"]

# chi = H_SO |S1>, sector blocks, filter polynomial (same choices as isc_filtered.filtered_isc)
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
psiA = A[:, 0]
p_true = float(np.sum(np.abs(psiA[good]) ** 2))
S_good = np.diag(np.where(good, -1.0, 1.0))
S_0 = np.eye(2 ** n)
S_0[0, 0] = -1.0
Q = -A @ S_0 @ A.conj().T @ S_good          # Grover operator; A|0> lives in a 2D invariant subspace
theta = math.asin(math.sqrt(p_true))
print(name, f"CAS({ncas},{ncas}) sector 2Sz={sz2_sel}: wires {n}, degree {deg}, p = {p_true:.3e}, theta = {theta:.4e}",
      flush=True)

# canonical QAE with t evaluation qubits, simulated with qml.QuantumPhaseEstimation
results = []
N_REP = 16
for t in range(4, 14):
    est = list(range(n, n + t))
    dev = qml.device("default.qubit", wires=n + t)

    @qml.qnode(dev)
    def qae():
        qml.QubitUnitary(A, wires=wires)
        qml.QuantumPhaseEstimation(qml.QubitUnitary(Q, wires=wires), estimation_wires=est)
        return qml.probs(wires=est)

    probs = np.asarray(qae())
    probs = np.clip(probs, 0, None) / np.sum(np.clip(probs, 0, None))
    y = np.arange(2 ** t)
    p_hat = np.sin(np.pi * y / 2 ** t) ** 2
    M = 2 ** t - 1                       # Grover calls in the QPE register (sum of 2^j)
    calls_run = 2 * M + 1               # calls to A or A^dagger in one QAE run
    calls = N_REP * calls_run           # median of N_REP runs
    rs = np.random.default_rng(1)
    draws = rs.choice(p_hat, size=(4000, N_REP), p=probs)
    med = np.median(draws, axis=1)
    err_ae = float(np.median(np.abs(med - p_true)))
    q90 = float(np.quantile(np.abs(med - p_true), 0.9))
    # direct sampling with the same number of calls: exact binomial statistics of the sample mean
    from scipy.stats import binom
    k = np.arange(0, min(calls, 20000) + 1)
    pk = binom.pmf(k, calls, p_true)
    errs = np.abs(k / calls - p_true)
    order = np.argsort(errs)
    cdf = np.cumsum(pk[order])
    err_samp = float(errs[order][np.searchsorted(cdf, 0.5)])
    q90_s = float(errs[order][min(np.searchsorted(cdf, 0.9), len(order) - 1)])
    results.append({"t": t, "calls": int(calls), "median_abs_err_ae": err_ae, "q90_abs_err_ae": q90,
                    "median_abs_err_sampling": err_samp, "q90_abs_err_sampling": q90_s,
                    "bound_ae": 2 * math.pi * math.sqrt(p_true) / M + math.pi ** 2 / M ** 2})
    print(f"t={t:2d} calls={calls:7d}  AE median|err|={err_ae:.2e} (90%: {q90:.2e})  "
          f"sampling median|err|={err_samp:.2e} (90%: {q90_s:.2e})", flush=True)

out = {"name": name, "ncas": ncas, "sz2": sz2_sel, "degree": deg, "p": p_true, "results": results}
os.makedirs(os.path.join(ROOT, "paper", "quantum", "results"), exist_ok=True)
json.dump(out, open(os.path.join(ROOT, "paper", "quantum", "results", "ae_demo.json"), "w"), indent=1)

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
fig, ax = plt.subplots(figsize=(3.4, 3.0))
x = [r["calls"] for r in results]
ax.loglog(x, [max(r["q90_abs_err_sampling"], 1e-30) / p_true for r in results], "s-", color="#c05621", ms=3,
          label="direct sampling")
ax.loglog(x, [max(r["q90_abs_err_ae"], 1e-30) / p_true for r in results], "o-", color="#2b6cb0", ms=3,
          label="amplitude estimation (median of 16)")
ax.set_xlabel("calls to the filter circuit")
ax.set_ylabel(r"90th-percentile relative error in $k_\mathrm{ISC}$")
ax.set_title(f"{name.replace('_', ' ')}, CAS({ncas},{ncas}), $p$ = {p_true:.1e}", fontsize=7)
ax.legend(frameon=False, fontsize=7)
fig.tight_layout()
fig.savefig(os.path.join(ROOT, "paper", "figures", "ae_demo.pdf"))
