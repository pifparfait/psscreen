"""Decision-aware fault-tolerant cost of the gap-resolved ISC observable.

The golden-rule filter (psscreen/isc_filtered.py) returns k_ISC = K p, where
    K = (2 pi / hbar) ||H_SO|S1>||^2 rho_max   and   p = <chi|P(H)^2|chi>/<chi|chi>
is the success probability of one QSVT circuit A (S1 preparation -> H_SO block encoding -> filter).
A screening campaign does not need p itself but a decision, so the relevant cost is the number of calls to A
needed to take that decision with confidence 1 - delta:

  D1  rate estimate: k_ISC to relative statistical error r (default r = 0.2, well inside the factor-2 budget);
  D2  triplet-yield class: is Phi_T > 0.5, i.e. k_ISC > k_r + k_ic?  The required absolute error on p is the
      distance |p - p_th| between the dye and the threshold.

Two estimators are compared:
  * direct sampling of the success flag (what arXiv:2512.15889 uses, Hoeffding/Chernoff bounds), and
  * amplitude estimation (Brassard et al. 2002): canonical QAE with M Grover iterations has error
    <= 2 pi sqrt(p(1-p))/M + pi^2/M^2 with probability >= 8/pi^2; the median of n_rep runs boosts this to
    1 - delta with n_rep = ceil(ln(1/delta) / (2 (8/pi^2 - 1/2)^2)).

Cost of one call to A (Toffolis), for CAS(12,11) with THC qubitization (see paper/quantum/thc_pennylane_resources.py):
  C_A = C_prep(S1) + C_SO + d * C_walk,
  C_walk = 2.97e3 (THC rank 22, angles in 5 batches, 158 logical qubits),
  d = c (lambda/Delta) sqrt(ln(1/eps_poly)) with c fitted to the preprint's d = 4.7571 lambda'/Delta + 321 at
      eps = 0.01 (Gaussian-filter scaling), Delta = sqrt(2 lambda_s k_B T) the MLJ line width,
  C_prep = sum-of-Slaters state preparation with D = 1e4 determinants (preprint's formula),
  C_SO  = one-body H_SO by basis rotation (3 components x N(N-1)/2 Givens rotations x rotation bits).
The filter centre needs E_S1: one QPE to the per-state precision of the accuracy budget (2 mHa) is added once.

Electronic inputs per dye (channels S1->T_n: gaps and SOC; k_r, k_ic) come from the TDDFT/DeltaSCF pipeline
results. ||chi||^2 is the sum of |V_n|^2 over the computed triplets (6 per dye), a lower bound on the full value,
so p is an upper bound and the costs are optimistic; a column with ||chi||^2 x 10 shows the sensitivity.

Usage (from the repository root): python paper/quantum/decision_cost.py
Outputs: paper/quantum/results/decision_cost.{json,md}, paper/figures/decision_cost.pdf
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from psscreen.rates import CM1_TO_EV, HBAR_EVS, KB_EV, MLJParams, fcwd_mlj  # noqa: E402

HA_EV = 27.211386245988
P = MLJParams()

# --------------------------------------------------------------------------- cost model (CAS(12,11), THC)
N_ORB = 11
LAMBDA = 15.9                   # Ha, one-norm (BODIPY core CAS(12,11))
C_WALK = 2.97e3                 # Toffolis per qubitized walk step (THC, 5 angle batches)
QUBITS = 158                    # logical qubits of that encoding
BETH = 27                       # rotation bits (preprint rule)
D_SOS = 1e4                     # determinants in the S1 sum-of-Slaters preparation
DELTA_HA = math.sqrt(2 * P.lambda_s * KB_EV * P.T) / HA_EV
C_FIT = 4.7571 / math.sqrt(math.log(100.0))
EPS_E_HA = 2.0e-3               # per-state QPE precision from the accuracy budget
DELTA_CONF = 0.05               # failure probability per decision
XANADU = {"isc_toffoli_per_shot": 7.96e7, "isc_shots": 483, "tp_toffoli_per_shot": 2.72e7, "tp_shots": 133}


def c_prep(D=D_SOS):
    L = math.ceil(math.log2(D))
    return (2 * L - 2) * D + 2 ** (L + 1) + D


def c_so(n=N_ORB):
    return 3 * n * (n - 1) / 2 * BETH


def degree(eps_poly):
    eps_poly = min(eps_poly, 0.01)        # never coarser than the preprint fit
    return C_FIT * (LAMBDA / DELTA_HA) * math.sqrt(math.log(1 / eps_poly)) + 321.2051


def c_call(eps_poly):
    return c_prep() + c_so() + degree(eps_poly) * C_WALK


def c_qpe_s1():
    steps = math.pi * LAMBDA / (2 * EPS_E_HA)
    return steps * C_WALK + c_prep()


# --------------------------------------------------------------------------- estimators
N_REP = math.ceil(math.log(1 / DELTA_CONF) / (2 * (8 / math.pi ** 2 - 0.5) ** 2))


def ae_iterations(p, eps):
    """Smallest M with 2 pi sqrt(p(1-p))/M + pi^2/M^2 <= eps."""
    a = 2 * math.pi * math.sqrt(p * (1 - p))
    # eps M^2 - a M - pi^2 >= 0
    return math.ceil((a + math.sqrt(a * a + 4 * eps * math.pi ** 2)) / (2 * eps))


def sampling_shots(p, eps):
    """Bernstein/Chernoff: shots for |p_hat - p| <= eps with prob 1 - delta (Bernoulli variance p(1-p))."""
    L = math.log(2 / DELTA_CONF)
    return math.ceil(2 * p * (1 - p) * L / eps ** 2 + 2 * L / (3 * eps))


# --------------------------------------------------------------------------- dyes
def dye_inputs(d):
    ph = d["photophysics"]
    ch = ph["channels"]
    chi2_cm2 = sum(c["V_cm"] ** 2 for c in ch)
    grid = np.linspace(-1.5, 0.5, 2001)
    rho_max = max(fcwd_mlj(float(x), P) for x in grid)
    K = 2 * math.pi / HBAR_EVS * chi2_cm2 * CM1_TO_EV ** 2 * rho_max
    k = ph["k_isc"]
    k_th = ph["k_r"] + ph["k_ic"]
    return dict(chi2_cm2=chi2_cm2, K=K, k=k, p=k / K, k_th=k_th, p_th=k_th / K)


def analyse_one(name, d, chi_scale=1.0, r=0.2):
    x = dye_inputs(d)
    K = x["K"] * chi_scale
    p, p_th = x["k"] / K, x["k_th"] / K
    out = {"name": name, "k_isc": x["k"], "k_threshold": x["k_th"], "K": K, "p": p, "p_th": p_th}
    # D1: relative error r on p
    eps1 = r * p
    epsp1 = max(r * math.sqrt(p) / 4, 1e-12)
    out["D1"] = {"eps_abs": eps1, "degree": degree(epsp1), "C_call": c_call(epsp1),
                 "shots": sampling_shots(p, eps1), "ae_M": ae_iterations(p, eps1), "ae_rep": N_REP}
    # D2: classification against p_th
    eps2 = abs(p - p_th)
    epsp2 = max(eps2 / (4 * math.sqrt(max(p, 1e-300))), 1e-12)
    out["D2"] = {"eps_abs": eps2, "degree": degree(epsp2), "C_call": c_call(epsp2),
                 "shots": sampling_shots(p, eps2), "ae_M": ae_iterations(p, eps2), "ae_rep": N_REP}
    for D in ("D1", "D2"):
        e = out[D]
        e["toffoli_sampling"] = e["shots"] * e["C_call"] + c_qpe_s1()
        e["toffoli_ae"] = e["ae_rep"] * (2 * e["ae_M"] + 1) * e["C_call"] + c_qpe_s1()
        e["speedup"] = e["toffoli_sampling"] / e["toffoli_ae"]
    return out


def load_all():
    runs = ROOT / "validation" / "runs"
    cand = ROOT / "screening" / "runs" / "candidates_round1"
    out = []
    for n in ["BODIPY_core", "BDP", "I-BDP", "I2-BDP", "I3-BDP", "azaBDP_Ph4"]:
        out.append((n, json.loads((runs / f"{n}.json").read_text())))
    for n in ["BDP-0", "I-BDP-0", "I2-BDP-0"]:
        out.append((n, json.loads((runs / "set3" / f"{n}_c0.json").read_text())))
    for n in ["Cl2-BDP", "Br2-BDP", "I1-BDP", "azaBDP_Br2", "azaBDP_BrPh2", "azaBDP_IPh2"]:
        out.append((n, json.loads((cand / f"{n}.json").read_text())))
    return out


def fmt(x):
    return f"{x:.1e}"


def main():
    dyes = load_all()
    rows = [analyse_one(n, d) for n, d in dyes]
    rows10 = {n: analyse_one(n, d, chi_scale=10.0) for n, d in dyes}
    res = {"model": {"lambda_Ha": LAMBDA, "C_walk": C_WALK, "qubits": QUBITS, "Delta_eV": DELTA_HA * HA_EV,
                     "C_prep": c_prep(), "C_SO": c_so(), "C_QPE_S1": c_qpe_s1(), "n_rep": N_REP,
                     "delta": DELTA_CONF, "xanadu": XANADU},
           "rows": rows, "rows_chi_x10": list(rows10.values())}
    outd = ROOT / "paper" / "quantum" / "results"
    outd.mkdir(parents=True, exist_ok=True)
    (outd / "decision_cost.json").write_text(json.dumps(res, indent=1))

    L = ["# Decision-aware cost of the gap-resolved ISC observable (generated by decision_cost.py)", "",
         f"CAS(12,11), THC qubitization: {QUBITS} logical qubits, {C_WALK:.2e} Toffolis per walk step, "
         f"λ = {LAMBDA} Ha, filter width Δ = {DELTA_HA*HA_EV:.3f} eV. "
         f"S1 preparation {c_prep():.1e}, H_SO encoding {c_so():.1e}, one QPE for E(S1) at 2 mHa "
         f"{c_qpe_s1():.1e} Toffolis. Confidence 1 − δ = {1-DELTA_CONF:.2f} per decision "
         f"(amplitude estimation: median of {N_REP} runs).", "",
         f"Reference, arXiv:2512.15889 (Orig-BDP, N = 11): ISC proxy {XANADU['isc_toffoli_per_shot']:.2e} Toffolis "
         f"× {XANADU['isc_shots']} shots = {XANADU['isc_toffoli_per_shot']*XANADU['isc_shots']:.1e}; threshold "
         f"projection {XANADU['tp_toffoli_per_shot']:.2e} × {XANADU['tp_shots']} = "
         f"{XANADU['tp_toffoli_per_shot']*XANADU['tp_shots']:.1e} Toffolis.", "",
         "## D1: k_ISC to 20% statistical error (dyes with k_ISC ≥ 1e6 s⁻¹; for slower dyes only D2 matters)", "",
         "| Dye | k_ISC (s⁻¹) | p | filter degree | Toffolis per call | shots (sampling) | Toffolis, sampling | "
         "Grover calls (AE) | Toffolis, AE | sampling/AE | Toffolis AE, ‖χ‖²×10 |",
         "|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        if r["k_isc"] < 1e6:
            continue
        e = r["D1"]
        L.append(f"| {r['name']} | {r['k_isc']:.1e} | {r['p']:.1e} | {e['degree']:.0f} | {fmt(e['C_call'])} | "
                 f"{e['shots']:.1e} | {fmt(e['toffoli_sampling'])} | {e['ae_rep']*(2*e['ae_M']+1):.1e} | "
                 f"{fmt(e['toffoli_ae'])} | {e['speedup']:.0f} | {fmt(rows10[r['name']]['D1']['toffoli_ae'])} |")
    L += ["", "## D2: is Φ_T above 0.5 (k_ISC > k_r + k_ic)?", "",
          "| Dye | k_ISC / threshold | p | |p − p_th| | shots (sampling) | Toffolis, sampling | Toffolis, AE | "
          "best of the two |", "|---|---|---|---|---|---|---|---|"]
    for r in rows:
        e = r["D2"]
        L.append(f"| {r['name']} | {r['k_isc']/r['k_threshold']:.1e} | {r['p']:.1e} | {e['eps_abs']:.1e} | "
                 f"{e['shots']:.1e} | {fmt(e['toffoli_sampling'])} | {fmt(e['toffoli_ae'])} | "
                 f"{fmt(min(e['toffoli_sampling'], e['toffoli_ae']))} |")
    (outd / "decision_cost.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))
    fig(rows)


def fig(rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ps = np.logspace(-12, -1, 200)
    fig_, ax = plt.subplots(figsize=(3.4, 3.0))
    tot_s, tot_a = [], []
    for p in ps:
        eps = 0.2 * p
        cc = c_call(max(0.2 * math.sqrt(p) / 4, 1e-12))
        tot_s.append(sampling_shots(p, eps) * cc + c_qpe_s1())
        tot_a.append(N_REP * (2 * ae_iterations(p, eps) + 1) * cc + c_qpe_s1())
    ax.loglog(ps, tot_s, color="#c05621", lw=1.2, label="direct sampling")
    ax.loglog(ps, tot_a, color="#2b6cb0", lw=1.2, label="amplitude estimation")
    for r in rows:
        if r["k_isc"] < 1e5:
            continue
        ax.plot(r["p"], r["D1"]["toffoli_ae"], "o", color="#2b6cb0", ms=3)
        ax.plot(r["p"], r["D1"]["toffoli_sampling"], "s", color="#c05621", ms=3)
    x = XANADU["isc_toffoli_per_shot"] * XANADU["isc_shots"]
    ax.axhline(x, color="0.5", ls=":", lw=0.8)
    ax.text(2e-12, x * 1.6, "gap-free proxy (preprint)", fontsize=6, color="0.4")
    ax.set_xlabel(r"success probability $p = k_\mathrm{ISC}/K$")
    ax.set_ylabel("Toffolis for $k_\\mathrm{ISC}$ to 20%")
    ax.legend(frameon=False, fontsize=7, loc="upper right")
    fig_.tight_layout()
    fig_.savefig(ROOT / "paper" / "figures" / "decision_cost.pdf")


if __name__ == "__main__":
    main()
