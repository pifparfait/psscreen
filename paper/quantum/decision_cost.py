"""Decision-aware fault-tolerant cost of the gap-resolved ISC observable.

The golden-rule filter (psscreen/isc_filtered.py) estimates k_ISC = (2 pi/hbar) <chi|rho(H - E_S1)|chi> with
chi = H_SO|S1>. On a fault-tolerant machine chi is produced by a block encoding of H_SO with normalization alpha_SO
(paper/quantum/alpha_so.py), so one call of the circuit A succeeds with probability
    p = k_ISC / K',   K' = (2 pi/hbar) alpha_SO^2 rho_max.
K' is known exactly from the Hamiltonian; no estimate of ||chi||^2 is needed. If K' < k_r + k_ic the triplet-yield
class is decided without any call (k_ISC <= K').

A screening campaign needs a decision per dye, taken with confidence 1 - delta:
  D1  rate estimate: k_ISC to relative statistical error r (default 0.2: the statistical part uses about a quarter
      of the factor-2 budget in ln k, leaving the rest to the systematic gap error);
  D2  triplet-yield class: is Phi_T > 0.5, i.e. k_ISC > k_r + k_ic?

Estimators (all with confidence 1 - delta = 0.95 per decision):
  * direct sampling of the success flag. D1: Dagum-Karp-Luby-Ross stopping rule, which needs no prior knowledge of
    p (expected calls Y1/p, Y1 = 1 + 4(e-2)(1+r) ln(2/delta)/r^2). D2: Wald's sequential probability ratio test
    against p_th, expected calls ~ ln(1/delta)/KL(p || p_th);
  * maximum-likelihood amplitude estimation (MLAE, no QPE and no controlled Grover operator): calls = c/(r sqrt(p))
    with the empirical 95th-percentile constant c from paper/quantum/ae_demo.py (simulated on the actual filter
    circuit, valid for r >= 0.05). For D2 the required angle precision is |asin sqrt(p_th) - asin sqrt(p)|;
  * canonical QAE (Brassard et al. 2002, Thm 12) as a rigorous check: M rounded up to a power of two, median of 16
    runs, 2(M-1)+1 calls per run; the overhead of controlling the walk operator is neglected.

Cost of one call (Toffolis), CAS(12,11) with THC qubitization (paper/quantum/thc_pennylane_resources.py):
  C_A = C_S1 + C_SO + d * C_walk
  C_walk = 2.97e3 (158 logical qubits, angles in 5 batches);
  d = filter degree, fitted to our own Chebyshev filters (paper/quantum/degree_check.py):
      d = s(eps) lambda'/Delta, s(eps) linear in ln(1/eps), lambda' = lambda + |E0 - c0| (shift of the
      block-encoded operator relative to its identity coefficient, 3.46 Ha for the BODIPY core), Delta the MLJ width;
      polynomial precision eps = r sqrt(p)/4 so that the induced error in p stays below r p / 2;
  C_S1 = coherent S1 preparation: sum-of-Slaters (D = 1e4 determinants) followed by an eigenstate filter of degree
      (lambda'/Delta_21) ln(1/eps_S) (Delta_21 = 0.5 eV, eps_S = 1e-3), boosted by one round of amplitude
      amplification (3 applications) for an assumed initial overlap gamma^2 = 0.5;
  C_SO = one-body H_SO by a single basis rotation (3 x N(N-1)/2 Givens rotations x rotation bits, conservative).
One QPE for E(S1) at the per-state precision of the accuracy budget (2 mHa) is added once per dye; its error enters
as a systematic gap error, not as statistics.

Electronic inputs per dye (k_ISC, k_r, k_ic) come from the TDDFT/DeltaSCF pipeline; alpha_SO from the CAS(12,11)
SOC integrals.

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
RES = ROOT / "paper" / "quantum" / "results"

# --------------------------------------------------------------------------- cost model (CAS(12,11), THC)
N_ORB = 11
LAMBDA = 15.9                   # Ha, one-norm (BODIPY core CAS(12,11))
_dc = json.loads((RES / "degree_check.json").read_text())
SHIFT = _dc["bodipy_core_cas12_11"]["abs_E0_minus_c0_Ha"]
LAMBDA_P = LAMBDA + SHIFT       # normalization of (H - E_c)
C_WALK = 2.97e3                 # Toffolis per qubitized walk step (THC, 5 angle batches)
QUBITS = 158                    # logical qubits of that encoding
BETH = 27                       # rotation bits (preprint rule)
D_SOS = 1e4                     # determinants in the S1 sum-of-Slaters preparation
GAMMA2 = 0.5                    # assumed overlap of the sum-of-Slaters state with S1
DELTA21_HA = 0.5 / HA_EV        # S1-S2 gap used for the S1 eigenstate filter
EPS_S = 1e-3
DELTA_HA = math.sqrt(2 * P.lambda_s * KB_EV * P.T) / HA_EV
EPS_E_HA = 2.0e-3               # per-state QPE precision from the accuracy budget
DELTA_CONF = 0.05               # failure probability per decision
R_D1 = 0.2
XANADU = {"isc_toffoli_per_shot": 7.96e7, "isc_shots": 483, "tp_toffoli_per_shot": 2.72e7, "tp_shots": 133}

# degree slope s(eps) = a ln(1/eps) + b from our own filters
_le = np.array([math.log(1 / float(k)) for k in _dc["fits"]])
_sl = np.array([v["slope"] for v in _dc["fits"].values()])
DEG_A, DEG_B = np.polyfit(_le, _sl, 1)
_ae = json.loads((RES / "ae_demo.json").read_text())
C_MLAE = _ae["c_mlae_q95"]      # eps_rel ~ C_MLAE / (sqrt(p) N_calls) at the confidence level of ae_demo


def c_prep(D=D_SOS):
    L = math.ceil(math.log2(D))
    return (2 * L - 2) * D + 2 ** (L + 1) + D


def amp_rounds(g2=GAMMA2):
    return max(0, math.ceil(math.pi / (4 * math.asin(math.sqrt(g2))) - 0.5))


def c_s1():
    d_s1 = LAMBDA_P / DELTA21_HA * math.log(1 / EPS_S)
    return (2 * amp_rounds() + 1) * (c_prep() + d_s1 * C_WALK)


def c_so(n=N_ORB):
    return 3 * n * (n - 1) / 2 * BETH


def degree(eps_poly):
    eps_poly = min(eps_poly, 0.01)
    return (DEG_A * math.log(1 / eps_poly) + DEG_B) * LAMBDA_P / DELTA_HA


def c_call(eps_poly):
    return c_s1() + c_so() + degree(eps_poly) * C_WALK


def c_qpe_s1():
    steps = math.pi * LAMBDA_P / (2 * EPS_E_HA)
    return steps * C_WALK + c_s1()


# --------------------------------------------------------------------------- estimators
N_REP = math.ceil(math.log(1 / DELTA_CONF) / (2 * (8 / math.pi ** 2 - 0.5) ** 2))


def canonical_calls(p, eps):
    """Canonical QAE: smallest power-of-two M with 2 pi sqrt(p(1-p))/M + pi^2/M^2 <= eps; median of N_REP runs."""
    a = 2 * math.pi * math.sqrt(p * (1 - p))
    M = (a + math.sqrt(a * a + 4 * eps * math.pi ** 2)) / (2 * eps)
    M = 2 ** math.ceil(math.log2(max(M, 1)))
    return N_REP * (2 * (M - 1) + 1)


def mlae_calls(p, r=None, p_th=None):
    """Grid-free amplitude estimation (MLAE) with the empirical constant c of eps_rel = c/(sqrt(p) N), i.e. an
    angle precision delta_theta = c/(2N). D1: delta_theta = r theta / 2. D2: delta_theta = |theta_th - theta|, which
    also covers p << p_th (certifying a small p needs resolution sqrt(p_th), not p_th)."""
    th = math.asin(math.sqrt(p))
    dth = r * th / 2 if p_th is None else abs(math.asin(math.sqrt(p_th)) - th)
    return max(C_MLAE / (2 * dth), 100.0)


def dklr_calls(p, r):
    """Dagum-Karp-Luby-Ross stopping rule for relative error r (expected number of Bernoulli trials)."""
    y1 = 1 + 4 * (math.e - 2) * (1 + r) * math.log(2 / DELTA_CONF) / r ** 2
    return y1 / p


def kl(a, b):
    return a * math.log(a / b) + (1 - a) * math.log((1 - a) / (1 - b))


def sprt_calls(p, p_th):
    """Wald SPRT, H0: p = p_th vs the true p, error delta on both sides; expected number of trials."""
    return max(math.log((1 - DELTA_CONF) / DELTA_CONF) / kl(p, p_th), 1.0)


# --------------------------------------------------------------------------- dyes
_alpha = json.loads((RES / "alpha_so.json").read_text())
GRID = np.linspace(-1.5, 0.5, 2001)
RHO_MAX = max(fcwd_mlj(float(x), P) for x in GRID)


def dye_inputs(name, d):
    ph = d["photophysics"]
    chi2_cm2 = sum(c["V_cm"] ** 2 for c in ph["channels"])
    a = _alpha[name]["alpha_so_cm"]
    Kp = 2 * math.pi / HBAR_EVS * (a * CM1_TO_EV) ** 2 * RHO_MAX
    k, k_th = ph["k_isc"], ph["k_r"] + ph["k_ic"]
    return dict(alpha_so_cm=a, chi_tddft_cm=math.sqrt(chi2_cm2), K=Kp, k=k, p=k / Kp, k_th=k_th, p_th=k_th / Kp)


def analyse_one(name, d, r=R_D1):
    x = dye_inputs(name, d)
    p, p_th = x["p"], x["p_th"]
    out = {"name": name, "k_isc": x["k"], "k_threshold": x["k_th"], "alpha_so_cm": x["alpha_so_cm"],
           "chi_tddft_cm": x["chi_tddft_cm"], "K": x["K"], "p": p, "p_th": p_th,
           "consistent": bool(p < 1)}
    q = qpe = c_qpe_s1()
    # D1
    if p < 1:
        epsp = r * math.sqrt(p) / 4
        cc = c_call(epsp)
        e = {"degree": degree(epsp), "C_call": cc, "calls_sampling": dklr_calls(p, r),
             "calls_mlae": mlae_calls(p, r=r), "calls_canonical": canonical_calls(p, r * p)}
        for k_ in ("sampling", "mlae", "canonical"):
            e["toffoli_" + k_] = e["calls_" + k_] * cc + qpe
        e["ratio_to_per_call"] = e["toffoli_mlae"] / cc
        out["D1"] = e
    # D2
    if p_th >= 1:
        out["D2"] = {"free": True, "toffoli_best": 0.0}
    elif p < 1:
        marg = abs(p - p_th) / p
        epsp = max(abs(p - p_th) / (4 * math.sqrt(p)), 1e-15)
        cc = c_call(epsp)
        e = {"free": False, "margin_rel": marg, "degree": degree(epsp), "C_call": cc,
             "calls_sampling": sprt_calls(p, p_th), "calls_mlae": mlae_calls(p, p_th=p_th),
             "calls_canonical": canonical_calls(p, abs(p - p_th))}
        for k_ in ("sampling", "mlae", "canonical"):
            e["toffoli_" + k_] = e["calls_" + k_] * cc + qpe
        e["toffoli_best"] = min(e["toffoli_sampling"], e["toffoli_mlae"])
        e["ratio_to_per_call"] = e["toffoli_best"] / cc
        out["D2"] = e
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
    model = {"lambda_Ha": LAMBDA, "lambda_prime_Ha": LAMBDA_P, "C_walk": C_WALK, "qubits": QUBITS,
             "Delta_eV": DELTA_HA * HA_EV, "degree_slope": [float(DEG_A), float(DEG_B)], "C_S1": c_s1(), "C_SO": c_so(),
             "C_QPE_S1": c_qpe_s1(), "n_rep_canonical": N_REP, "c_mlae": C_MLAE, "delta": DELTA_CONF,
             "gamma2": GAMMA2, "r_D1": R_D1, "xanadu": XANADU}
    (RES / "decision_cost.json").write_text(json.dumps({"model": model, "rows": rows}, indent=1))

    L = ["# Decision-aware cost of the gap-resolved ISC observable (generated by decision_cost.py)", "",
         f"CAS(12,11), THC qubitization: {QUBITS} logical qubits, {C_WALK:.2e} Toffolis per walk step, "
         f"λ' = λ + |E0 − c0| = {LAMBDA} + {SHIFT:.2f} = {LAMBDA_P:.1f} Ha, filter width Δ = {DELTA_HA*HA_EV:.3f} eV, "
         f"filter degree d = ({DEG_A:.3f} ln(1/ε) + {DEG_B:.3f}) λ'/Δ (fitted to our filters). "
         f"Coherent S1 preparation (sum of Slaters + eigenstate filter, γ² = {GAMMA2}, {2*amp_rounds()+1} "
         f"applications) {c_s1():.1e}, H_SO basis rotation {c_so():.1e}, one QPE for E(S1) at 2 mHa "
         f"{c_qpe_s1():.1e} Toffolis. Success probability p = k_ISC/K′, K′ = (2π/ħ) α_SO² ρ_max. "
         f"Confidence 0.95 per decision. MLAE constant c = {C_MLAE:.1f} (95th percentile, from ae_demo.py).", "",
         f"Reference, arXiv:2512.15889 (Orig-BDP, N = 11): ISC proxy {XANADU['isc_toffoli_per_shot']:.2e} Toffolis "
         f"× {XANADU['isc_shots']} shots = {XANADU['isc_toffoli_per_shot']*XANADU['isc_shots']:.1e} at fixed additive "
         "precision (not a decision; the proxy gives no absolute rate).", "",
         "## Spin-orbit normalization", "",
         "| Dye | α_SO, CAS(12,11) (cm⁻¹) | ‖H_SO S1‖, TDDFT (cm⁻¹) | k_ISC (s⁻¹) | K′ (s⁻¹) | p = k_ISC/K′ |",
         "|---|---|---|---|---|---|"]
    for r in rows:
        flag = "" if r["consistent"] else " (> 1: rate not reachable in this active space)"
        L.append(f"| {r['name']} | {r['alpha_so_cm']:.0f} | {r['chi_tddft_cm']:.0f} | {r['k_isc']:.1e} | "
                 f"{r['K']:.1e} | {r['p']:.1e}{flag} |")
    L += ["", f"## D1: k_ISC to {int(R_D1*100)}% statistical error (dyes with k_ISC ≥ 1e6 s⁻¹; slower dyes have "
          "Φ_T far below 0.5, so only D2 matters)", "",
          "| Dye | p | degree | Toffolis per call | calls, sampling (DKLR) | calls, MLAE | calls, canonical QAE | "
          "Toffolis, sampling | Toffolis, MLAE | sampling/MLAE | total/per call |",
          "|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        if r["k_isc"] < 1e6 or "D1" not in r:
            continue
        e = r["D1"]
        L.append(f"| {r['name']} | {r['p']:.1e} | {e['degree']:.0f} | {fmt(e['C_call'])} | "
                 f"{e['calls_sampling']:.1e} | {e['calls_mlae']:.1e} | {e['calls_canonical']:.1e} | "
                 f"{fmt(e['toffoli_sampling'])} | {fmt(e['toffoli_mlae'])} | "
                 f"{e['toffoli_sampling']/e['toffoli_mlae']:.0f} | {e['ratio_to_per_call']:.1e} |")
    L += ["", "## D2: is Φ_T above 0.5 (k_ISC > k_r + k_ic)?", "",
          "| Dye | k_ISC / threshold | p | p_th | calls, sampling (SPRT) | calls, MLAE | Toffolis, sampling | "
          "Toffolis, MLAE | best | total/per call |", "|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        e = r.get("D2")
        if e is None:
            continue
        if e["free"]:
            L.append(f"| {r['name']} | {r['k_isc']/r['k_threshold']:.1e} | {r['p']:.1e} | {r['p_th']:.1e} | "
                     "0 (K′ < threshold) | 0 | 0 | 0 | 0 | 0 |")
            continue
        L.append(f"| {r['name']} | {r['k_isc']/r['k_threshold']:.1e} | {r['p']:.1e} | {r['p_th']:.1e} | "
                 f"{e['calls_sampling']:.1e} | {e['calls_mlae']:.1e} | {fmt(e['toffoli_sampling'])} | "
                 f"{fmt(e['toffoli_mlae'])} | {fmt(e['toffoli_best'])} | {e['ratio_to_per_call']:.1e} |")
    L += ["", f"δ = {DELTA_CONF} per decision. For a campaign of N dyes with a family-wise error δ, replace δ by δ/N: "
          "the cost grows by ln(2N/δ)/ln(2/δ), e.g. ×2.2 for N = 100. With r = 0.4 instead of 0.2 the D1 MLAE "
          "costs halve."]
    (RES / "decision_cost.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))
    fig(rows)


def fig(rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ps = np.logspace(-12, -0.5, 200)
    fig_, ax = plt.subplots(figsize=(3.4, 3.0))
    tot_s, tot_a = [], []
    for p in ps:
        cc = c_call(R_D1 * math.sqrt(p) / 4)
        tot_s.append(dklr_calls(p, R_D1) * cc + c_qpe_s1())
        tot_a.append(mlae_calls(p, r=R_D1) * cc + c_qpe_s1())
    ax.loglog(ps, tot_s, color="#c05621", lw=1.2, label="direct sampling")
    ax.loglog(ps, tot_a, color="#2f855a", lw=1.2, label="amplitude estimation (MLAE)")
    for r in rows:
        if r["k_isc"] < 1e6 or "D1" not in r:
            continue
        ax.plot(r["p"], r["D1"]["toffoli_mlae"], "o", color="#2f855a", ms=3)
        ax.plot(r["p"], r["D1"]["toffoli_sampling"], "s", color="#c05621", ms=3)
    x = XANADU["isc_toffoli_per_shot"] * XANADU["isc_shots"]
    ax.axhline(x, color="0.5", ls=":", lw=0.8)
    ax.text(2e-12, x * 1.6, "proxy, fixed additive precision (not a decision)", fontsize=5.5, color="0.4")
    ax.set_xlabel(r"success probability $p = k_\mathrm{ISC}/K'$")
    ax.set_ylabel(f"Toffolis for $k_\\mathrm{{ISC}}$ to {int(R_D1*100)}%")
    ax.legend(frameon=False, fontsize=7, loc="upper right")
    fig_.tight_layout()
    fig_.savefig(ROOT / "paper" / "figures" / "decision_cost.pdf")


if __name__ == "__main__":
    main()
