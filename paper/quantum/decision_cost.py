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
  D2  triplet-yield class: is Phi_T > 0.5, i.e. k_ISC > k_r + k_ic? Tested with an indifference zone
      k_th/F < k_ISC < k_th F (F = 2, the factor-2 budget) and no knowledge of the margin: inside the zone the test
      may return either class (or 'within a factor of 2'). The first version, with the margin known in advance (an
      oracle margin), is kept in the json under D2.oracle.

Estimators (all with confidence 1 - delta = 0.95 per decision):
  * direct sampling of the success flag. D1: Dagum-Karp-Luby-Ross stopping rule, which needs no prior knowledge of
    p (expected calls Y1/p, Y1 = 1 + 4(e-2)(1+r) ln(2/delta)/r^2). D2: Wald's sequential probability ratio test
    between the zone edges p_th/F and p_th F, expected calls from Wald's approximations at the true p;
  * maximum-likelihood amplitude estimation (MLAE, no QPE and no controlled Grover operator): calls = c/(r sqrt(p))
    with the empirical 95th-percentile constant c from paper/quantum/ae_demo.py (simulated on the actual filter
    circuit, valid for r >= 0.05). For D2 the precision is halved stage by stage until the interval excludes one zone edge (final angle
    precision max(theta - theta0, theta1 - theta)/2, at most 2x the last stage in total, confidence delta/L per stage
    with the constant scaled by the normal quantile);
  * canonical QAE (Brassard et al. 2002, Thm 12) as a rigorous check: M rounded up to a power of two, median of 16
    runs, 2(M-1)+1 calls per run; the overhead of controlling the walk operator is neglected.

Cost of one call (Toffolis), CAS(12,11) with THC qubitization (paper/quantum/thc_pennylane_resources.py):
  C_A = C_S1 + C_SO + d * C_walk
  C_walk = 3.25e3 (THC rank 33, 193 logical qubits, angles in 5 batches);
  d = filter degree of our own Chebyshev filters: d = s(eps) lambda'/Delta, with s(eps) computed directly for
      eps = 1e-2 ... 1e-7 (paper/quantum/degree_direct.py; linear in alpha/Delta, degree_check.py) and interpolated
      in ln(1/eps), lambda = 17.3 Ha (regularized THC fit), lambda' = lambda + |E0 - c0| (shift of the
      block-encoded operator relative to its identity coefficient, 3.46 Ha for the BODIPY core), Delta the MLJ width;
      polynomial precision eps = r sqrt(p)/4 so that the induced error in p stays below r p / 2;
  C_S1 = coherent S1 preparation: a spin-adapted sum of D = 1e4 configuration state functions (so that triplet
      components, which a spin-free eigenstate filter cannot reject, are absent by construction) followed by an
      eigenstate filter of degree (lambda'/Delta_21) ln(1/eps_S) (Delta_21 = 0.5 eV, eps_S = 1e-3). The filter flag
      succeeds with probability gamma^2 (assumed 0.5); it is folded into the amplitude-estimation good state, so the
      success probability seen by the estimators is gamma^2 p and the call count scales as 1/gamma;
  C_SO = one-body H_SO by a single basis rotation (3 x N(N-1)/2 Givens rotations x rotation bits, conservative).
One QPE for E(S1) at 0.1 mHa (0.0027 eV, 1/30 of the 0.08 eV gap budget and 1/6 of the 0.016 eV near-resonant
tolerance) is added once per dye; its error enters as a systematic gap error, not as statistics.
A channel-resolved baseline is also costed: QPE for each triplet channel that carries >= 1% of k_ISC, plus amplitude
estimation of V_n^2 through a T_n eigenstate filter (the same degree model with the line width replaced by the
distance to the nearest other triplet), with the golden rule evaluated classically.

Electronic inputs per dye (k_ISC, k_r, k_ic) come from the TDDFT/DeltaSCF pipeline; alpha_SO from the CAS(12,11)
SOC integrals.

Usage (from the repository root): python paper/quantum/decision_cost.py
Outputs: paper/quantum/results/decision_cost.{json,md}, paper/figures/decision_cost.pdf
"""
from __future__ import annotations

import json
import math
import os
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
# Ha, one-norm of the regularized rank-33 THC fit of the BODIPY core CAS(12,11) (paper/quantum/thc_regularized.py;
# double factorization gives 15.9). Set PSSCREEN_LAMBDA to test another value, e.g. 37.5 (unregularized THC fit,
# paper/quantum/thc_error.py); outputs then get a _lambda<value> suffix.
LAMBDA = float(os.environ.get("PSSCREEN_LAMBDA", 17.3))
SUFFIX = "" if "PSSCREEN_LAMBDA" not in os.environ else f"_lambda{LAMBDA:g}"
_dc = json.loads((RES / "degree_check.json").read_text())
SHIFT = _dc["bodipy_core_cas12_11"]["abs_E0_minus_c0_Ha"]
LAMBDA_P = LAMBDA + SHIFT       # normalization of (H - E_c)
C_WALK = 3.25e3                 # Toffolis per qubitized walk step (THC rank 3N = 33, 5 angle batches; rank 2N fails the gap budget, paper/quantum/thc_error.py)
QUBITS = 193                    # logical qubits of that encoding
BETH = 27                       # rotation bits (preprint rule)
D_SOS = 1e4                     # determinants in the S1 sum-of-Slaters preparation
GAMMA2 = 0.5                    # assumed overlap of the sum-of-Slaters state with S1
DELTA21_HA = 0.5 / HA_EV        # S1-S2 gap used for the S1 eigenstate filter (spin-adapted trial state)
EPS_S = 1e-3
DELTA_HA = math.sqrt(2 * P.lambda_s * KB_EV * P.T) / HA_EV
EPS_E_HA = 0.1e-3               # QPE precision for the filter centre: 0.0027 eV, 1/30 of the 0.08 eV budget and 1/6 of
                                # the 0.016 eV near-resonant tolerance
F_IND = 2.0                     # D2 indifference zone: k within a factor F_IND of the threshold is not decided
TOFFOLI_RATE = 1e6              # Toffolis per second, for wall-clock estimates
DELTA_CONF = 0.05               # failure probability per decision
R_D1 = 0.2
XANADU = {"isc_toffoli_per_shot": 7.96e7, "isc_shots": 483, "tp_toffoli_per_shot": 2.72e7, "tp_shots": 133}

# degree slope s(eps) = d / (alpha/Delta), computed directly on a grid eps = 1e-2 ... 1e-7 (degree_direct.py),
# interpolated in ln(1/eps); below 1e-7 linear extrapolation from the last two points
_dd = json.loads((RES / "degree_direct.json").read_text())
_rk = "400" if len(_dd["ratios"].get("400", {}).get("slope", [])) == len(_dd["eps"]) else "200"
_LE = np.array([math.log(1 / e) for e in _dd["eps"]])
_SL = np.array(_dd["ratios"][_rk]["slope"])
DEG_A, DEG_B = np.polyfit(_LE, _SL, 1)    # reported only (linear summary of the grid)


def slope(eps_poly):
    le = math.log(1 / min(eps_poly, 0.01))
    if le <= _LE[-1]:
        return float(np.interp(le, _LE, _SL))
    return float(_SL[-1] + (_SL[-1] - _SL[-2]) / (_LE[-1] - _LE[-2]) * (le - _LE[-1]))


_ae = json.loads((RES / "ae_demo.json").read_text())
C_MLAE = _ae["c_mlae_q95"]      # eps_rel ~ C_MLAE / (sqrt(p) N_calls) at the confidence level of ae_demo


def c_prep(D=D_SOS):
    L = math.ceil(math.log2(D))
    return (2 * L - 2) * D + 2 ** (L + 1) + D


def c_s1(gap_ha=DELTA21_HA):
    d_s1 = LAMBDA_P / gap_ha * math.log(1 / EPS_S)
    return c_prep() + d_s1 * C_WALK


def c_so(n=N_ORB):
    return 3 * n * (n - 1) / 2 * BETH


def degree(eps_poly, width_ha=DELTA_HA):
    return slope(eps_poly) * LAMBDA_P / width_ha


def c_call(eps_poly):
    return c_s1() + c_so() + degree(eps_poly) * C_WALK


def c_qpe_s1(eps_ha=EPS_E_HA):
    steps = math.pi * LAMBDA_P / (2 * eps_ha)
    return (steps * C_WALK + c_s1()) / GAMMA2      # repeated until the S1 flag succeeds


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
    """Wald SPRT, H0: p = p_th vs the true p, error delta on both sides; expected number of trials (oracle margin)."""
    return max(math.log((1 - DELTA_CONF) / DELTA_CONF) / kl(p, p_th), 1.0)


def sprt_zone_calls(p, p0, p1, delta=DELTA_CONF):
    """Wald SPRT of H0: p <= p0 against H1: p >= p1 (indifference zone p0 < p < p1), errors delta, evaluated at the
    true p with Wald's approximations for the operating characteristic L(p) and the average sample number. Needs no
    knowledge of the margin."""
    lnA, lnB = math.log((1 - delta) / delta), math.log(delta / (1 - delta))
    zs, zf = math.log(p1 / p0), math.log((1 - p1) / (1 - p0))
    drift = p * zs + (1 - p) * zf
    var = p * zs ** 2 + (1 - p) * zf ** 2 - drift ** 2
    if abs(drift) < 1e-12 * var ** 0.5:
        return -lnA * lnB / var
    # h != 0 with p (p1/p0)^h + (1-p) ((1-p1)/(1-p0))^h = 1, on the side opposite to the drift
    from scipy.optimize import brentq
    g = lambda h: p * math.expm1(h * zs) + (1 - p) * math.expm1(h * zf)  # noqa: E731
    sgn = 1.0 if drift < 0 else -1.0
    H = 1.0
    while g(sgn * H) < 0:
        H *= 2
    h = brentq(g, sgn * 1e-9 * H, sgn * H) if sgn > 0 else brentq(g, sgn * H, sgn * 1e-9 * H)
    Ah, Bh = math.expm1(h * lnA), math.expm1(h * lnB)
    L = Ah / (Ah - Bh)          # probability of accepting H0
    return max((L * lnB + (1 - L) * lnA) / drift, 1.0)


def z_two_sided(delta):
    from scipy.stats import norm
    return float(norm.isf(delta / 2))


def adaptive_stages(theta, dth):
    """Number of precision-halving stages from dth_0 = theta (or pi/4) down to dth."""
    return max(1, math.ceil(math.log2(max(min(theta, math.pi / 4), dth) / dth)) + 1)


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
    return dict(alpha_so_cm=a, chi_tddft_cm=math.sqrt(chi2_cm2), K=Kp, k=k, p=k / Kp, k_th=k_th, p_th=k_th / Kp,
                v_eff_cm=a * math.sqrt(k / Kp))


def canonical_depth(p, eps):
    """Largest Grover power of one canonical-QAE run (M - 1)."""
    a = 2 * math.pi * math.sqrt(p * (1 - p))
    M = (a + math.sqrt(a * a + 4 * eps * math.pi ** 2)) / (2 * eps)
    return 2 ** math.ceil(math.log2(max(M, 1))) - 1


def channel_route(d, alpha_cm, r=R_D1):
    """QPE per dominant channel + AE of V_n^2 through a T_n eigenstate filter; golden rule classical."""
    ch = d["photophysics"]["channels"]
    k = d["photophysics"]["k_isc"]
    dom = [c for c in ch if c["k"] >= 0.01 * k]
    tot, parts = 0.0, []
    for c in dom:
        others = [abs(c["dE_eV"] - o["dE_eV"]) for o in ch if o is not c]
        gap = max(min(others) if others else 0.5, 0.05) / HA_EV
        pn = GAMMA2 * (c["V_cm"] / alpha_cm) ** 2
        rn = r * k / c["k"] if len(dom) > 1 else r      # share of the relative error
        rn = min(rn, 0.5)
        eps = rn * math.sqrt(pn) / 4
        d_tn = degree(eps, gap)               # same degree model, with the line width replaced by the gap
        cc = c_s1() + c_so() + d_tn * C_WALK
        calls = mlae_calls(pn, r=rn)
        cost = calls * cc + c_qpe_s1()
        parts.append({"T": c["T"], "p_n": pn, "gap_eV": gap * HA_EV, "calls": calls, "C_call": cc, "toffoli": cost})
        tot += cost
    return {"toffoli": tot + c_qpe_s1(), "channels": parts}


def analyse_one(name, d, r=R_D1):
    x = dye_inputs(name, d)
    p_raw = x["p"]
    p, p_th = GAMMA2 * x["p"], GAMMA2 * x["p_th"]       # S1 filter flag folded into the good state
    out = {"name": name, "k_isc": x["k"], "k_threshold": x["k_th"], "alpha_so_cm": x["alpha_so_cm"],
           "chi_tddft_cm": x["chi_tddft_cm"], "v_eff_cm": x["v_eff_cm"], "K": x["K"], "p_raw": p_raw,
           "p": p, "p_th": p_th, "consistent": bool(p_raw < 1)}
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
        e["max_depth_toffoli"] = (2 * canonical_depth(p, r * p) + 1) * cc
        e["days_mlae"] = e["toffoli_mlae"] / TOFFOLI_RATE / 86400
        e["days_canonical"] = e["toffoli_canonical"] / TOFFOLI_RATE / 86400
        if x["k"] >= 1e6:
            e["channel_route"] = channel_route(d, x["alpha_so_cm"], r)
        out["D1"] = e
    # D2
    if x["p_th"] >= 1:
        out["D2"] = {"free": True, "toffoli_best": 0.0}
    elif p < 1:
        marg = abs(p - p_th) / p
        # (a) oracle margin, as in the first version: resolution |theta - theta_th|
        epsp = max(abs(p - p_th) / (4 * math.sqrt(p)), 1e-15)
        cc = c_call(epsp)
        e = {"free": False, "margin_rel": marg, "oracle": {
            "degree": degree(epsp), "C_call": cc, "calls_sampling": sprt_calls(p, p_th),
            "calls_mlae": mlae_calls(p, p_th=p_th), "calls_canonical": canonical_calls(p, abs(p - p_th))}}
        for k_ in ("sampling", "mlae", "canonical"):
            e["oracle"]["toffoli_" + k_] = e["oracle"]["calls_" + k_] * cc + qpe
        # (b) indifference zone k_th/F < k < k_th F, margin unknown. Sampling: SPRT between the zone edges.
        # Amplitude estimation: precision halved stage by stage until the confidence interval excludes one edge,
        # i.e. final angle precision dth = max(theta - theta0, theta1 - theta)/2; each stage at confidence
        # delta/L (L stages), the MLAE constant scaled by the two-sided normal quantile z(delta/L)/z(delta);
        # total calls <= 2x the last stage (geometric schedule).
        p0, p1 = p_th / F_IND, min(p_th * F_IND, 0.999)
        th, th0, th1 = (math.asin(math.sqrt(v)) for v in (p, p0, p1))
        dth = max(th - th0, th1 - th) / 2
        L = adaptive_stages(th, dth)
        zfac = z_two_sided(DELTA_CONF / L) / z_two_sided(DELTA_CONF)
        eps_z = max(dth / 2, 1e-15)                 # amplitude error of the polynomial below dth/2
        ccz = c_call(eps_z)
        mlae_z = 2 * max(C_MLAE * zfac / (2 * dth), 100.0)
        nrep_L = math.ceil(math.log(L / DELTA_CONF) / (2 * (8 / math.pi ** 2 - 0.5) ** 2))
        dp = math.sin(th + dth) ** 2 - p
        canon_z = 2 * canonical_calls(p, dp) * nrep_L / N_REP
        inside = bool(p0 < p < p1)
        e["zone"] = {"F": F_IND, "inside": inside, "stages": L, "z_factor": zfac, "degree": degree(eps_z),
                     "C_call": ccz, "calls_sampling": sprt_zone_calls(p, p0, p1), "calls_mlae": mlae_z,
                     "calls_canonical": canon_z}
        for k_ in ("sampling", "mlae", "canonical"):
            e["zone"]["toffoli_" + k_] = e["zone"]["calls_" + k_] * ccz + qpe
        # headline numbers: the indifference-zone test (no oracle margin)
        z = e["zone"]
        e.update({"degree": z["degree"], "C_call": ccz, "calls_sampling": z["calls_sampling"],
                  "calls_mlae": z["calls_mlae"], "calls_canonical": z["calls_canonical"],
                  "toffoli_sampling": z["toffoli_sampling"], "toffoli_mlae": z["toffoli_mlae"],
                  "toffoli_canonical": z["toffoli_canonical"], "inside_zone": inside})
        e["toffoli_best"] = min(e["toffoli_sampling"], e["toffoli_mlae"])
        e["toffoli_oracle_best"] = min(e["oracle"]["toffoli_sampling"], e["oracle"]["toffoli_mlae"])
        e["days_best"] = e["toffoli_best"] / TOFFOLI_RATE / 86400
        e["ratio_to_per_call"] = e["toffoli_best"] / ccz
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
    (RES / f"decision_cost{SUFFIX}.json").write_text(json.dumps({"model": model, "rows": rows}, indent=1))

    L = ["# Decision-aware cost of the gap-resolved ISC observable (generated by decision_cost.py)", "",
         f"CAS(12,11), THC qubitization: {QUBITS} logical qubits, {C_WALK:.2e} Toffolis per walk step, "
         f"λ' = λ + |E0 − c0| = {LAMBDA} + {SHIFT:.2f} = {LAMBDA_P:.1f} Ha, filter width Δ = {DELTA_HA*HA_EV:.3f} eV, "
         f"filter degree d = ({DEG_A:.3f} ln(1/ε) + {DEG_B:.3f}) λ'/Δ (fitted to our filters). "
         f"Coherent S1 preparation (spin-adapted sum of CSFs + eigenstate filter, flag success γ² = {GAMMA2} folded into "
         f"the good state) {c_s1():.1e}, H_SO basis rotation {c_so():.1e}, one QPE for E(S1) at "
         f"{EPS_E_HA*1e3:.1f} mHa {c_qpe_s1():.1e} Toffolis. Success probability γ² k_ISC/K′, "
         f"K′ = (2π/ħ) α_SO² ρ_max. "
         f"Confidence 0.95 per decision. MLAE constant c = {C_MLAE:.1f} (95th percentile, from ae_demo.py).", "",
         f"Reference, arXiv:2512.15889 (Orig-BDP, N = 11): ISC proxy {XANADU['isc_toffoli_per_shot']:.2e} Toffolis "
         f"× {XANADU['isc_shots']} shots = {XANADU['isc_toffoli_per_shot']*XANADU['isc_shots']:.1e} at fixed additive "
         "precision (not a decision; the proxy gives no absolute rate).", "",
         "## Spin-orbit normalization", "",
         "| Dye | α_SO, CAS(12,11) (cm⁻¹) | ‖H_SO S1‖, TDDFT (cm⁻¹) | V_eff = α_SO √(k/K′) (cm⁻¹) | k_ISC (s⁻¹) | "
         "K′ (s⁻¹) | k_ISC/K′ | p = γ² k_ISC/K′ |",
         "|---|---|---|---|---|---|---|---|"]
    for r in rows:
        flag = "" if r["consistent"] else " (> 1: rate not reachable in this active space)"
        L.append(f"| {r['name']} | {r['alpha_so_cm']:.0f} | {r['chi_tddft_cm']:.0f} | {r['v_eff_cm']:.2g} | "
                 f"{r['k_isc']:.1e} | {r['K']:.1e} | {r['p_raw']:.1e}{flag} | {r['p']:.1e} |")
    L += ["", f"## D1: k_ISC to {int(R_D1*100)}% statistical error (dyes with k_ISC ≥ 1e6 s⁻¹; slower dyes have "
          "Φ_T far below 0.5, so only D2 matters)", "",
          "| Dye | p | degree | Toffolis per call | calls, sampling (DKLR) | calls, MLAE | calls, canonical QAE | "
          "Toffolis, sampling | Toffolis, MLAE | Toffolis, canonical | sampling/MLAE | total/per call | "
          "deepest circuit | days at 1e6 Toffoli/s (MLAE / canonical) | channel-resolved route |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        if r["k_isc"] < 1e6 or "D1" not in r:
            continue
        e = r["D1"]
        L.append(f"| {r['name']} | {r['p']:.1e} | {e['degree']:.0f} | {fmt(e['C_call'])} | "
                 f"{e['calls_sampling']:.1e} | {e['calls_mlae']:.1e} | {e['calls_canonical']:.1e} | "
                 f"{fmt(e['toffoli_sampling'])} | {fmt(e['toffoli_mlae'])} | {fmt(e['toffoli_canonical'])} | "
                 f"{e['toffoli_sampling']/e['toffoli_mlae']:.0f} | {e['ratio_to_per_call']:.1e} | "
                 f"{fmt(e['max_depth_toffoli'])} | {e['days_mlae']:.0f} / {e['days_canonical']:.0f} | "
                 f"{fmt(e['channel_route']['toffoli'])} ({len(e['channel_route']['channels'])} ch.) |")
    L += ["", "## D2: is Φ_T above 0.5 (k_ISC > k_r + k_ic)?", "",
          "| Dye | k_ISC / threshold | p | p_th | calls, sampling (SPRT) | calls, MLAE | Toffolis, sampling | "
          "Toffolis, MLAE | Toffolis, canonical | best | inside zone | best, oracle margin | total/per call |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        e = r.get("D2")
        if e is None:
            continue
        if e["free"]:
            L.append(f"| {r['name']} | {r['k_isc']/r['k_threshold']:.1e} | {r['p']:.1e} | {r['p_th']:.1e} | "
                     "0 (K′ < threshold) | 0 | 0 | 0 | 0 | 0 | no | 0 | 0 |")
            continue
        L.append(f"| {r['name']} | {r['k_isc']/r['k_threshold']:.1e} | {r['p']:.1e} | {r['p_th']:.1e} | "
                 f"{e['calls_sampling']:.1e} | {e['calls_mlae']:.1e} | {fmt(e['toffoli_sampling'])} | "
                 f"{fmt(e['toffoli_mlae'])} | {fmt(e['toffoli_canonical'])} | {fmt(e['toffoli_best'])} | "
                 f"{'yes' if e['inside_zone'] else 'no'} | {fmt(e['toffoli_oracle_best'])} | {e['ratio_to_per_call']:.1e} |")
    nrep100 = math.ceil(math.log(100 / DELTA_CONF) / (2 * (8 / math.pi ** 2 - 0.5) ** 2))
    L += ["", f"δ = {DELTA_CONF} per decision. For a campaign of N dyes with a family-wise error δ, replace δ by δ/N. "
          f"For canonical QAE (median of runs) this raises the number of runs from {N_REP} to {nrep100} for N = 100 "
          f"(×{nrep100/N_REP:.1f}); for MLAE the tail has not been characterized. D2 uses an indifference zone of a factor "
          f"{F_IND:g} around the threshold and no knowledge of the margin (dyes marked 'inside' may return either class); "
          "the oracle-margin costs of the first version are kept in the json (D2.oracle). "
          "With r = 0.4 instead of 0.2 the D1 MLAE costs halve. Sensitivity: γ² = 0.1 / 0.9 multiplies the calls by "
          f"{math.sqrt(0.5/0.1):.2f} / {math.sqrt(0.5/0.9):.2f}; an S1 filter gap of 0.1 eV instead of 0.5 eV "
          f"(no spin-adapted trial state) raises C_S1 from {c_s1():.1e} to {c_s1(0.1/HA_EV):.1e}."]
    (RES / f"decision_cost{SUFFIX}.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))
    if not SUFFIX:
        fig(rows)
        latex_table(rows)


TABLE_DYES = ["BDP", "I-BDP", "I2-BDP", "I3-BDP", "I2-BDP-0", "Br2-BDP", "azaBDP_Br2", "azaBDP_IPh2"]
TEX_NAME = {"azaBDP_Br2": "aza-BDP-Br$_2$", "azaBDP_IPh2": "aza-BDP-IPh$_2$"}


def _tex(x, sig=2):
    if x == 0:
        return "0"
    e = int(math.floor(math.log10(abs(x))))
    m = x / 10 ** e
    if round(m, sig - 1) >= 10:
        m, e = m / 10, e + 1
    return f"${m:.{sig-1}f}\\times10^{{{e}}}$"


def _plain(x):
    if 0.01 <= x < 1:
        return f"{x:.2f}"
    if 1 <= x < 10:
        return f"{x:.1f}"
    if 10 <= x < 1000:
        return f"{x:.0f}"
    return _tex(x, 1)


def latex_table(rows):
    """Body rows of Table 2 of the paper (paper/quantum/results/decision_table.tex)."""
    by = {r["name"]: r for r in rows}
    out = []
    for n in TABLE_DYES:
        r = by[n]
        d1 = r.get("D1") if r["k_isc"] >= 1e6 else None
        e2 = r["D2"]
        d2 = "0" if e2.get("free") else _tex(e2["toffoli_best"]) + ("$^a$" if e2["inside_zone"] else "")
        cells = [TEX_NAME.get(n, n), f"{r['alpha_so_cm']:.0f}", _plain(r["v_eff_cm"]),
                 _plain(r["k_isc"] / r["k_threshold"]), _tex(r["p"]),
                 _tex(d1["toffoli_mlae"]) if d1 else "--", _tex(d1["toffoli_canonical"]) if d1 else "--",
                 _tex(d1["channel_route"]["toffoli"]) if d1 else "--", d2]
        out.append("    " + " & ".join(cells) + " \\\\")
    (RES / "decision_table.tex").write_text("% generated by paper/quantum/decision_cost.py\n" + "\n".join(out) + "\n")


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
