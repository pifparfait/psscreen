"""Multimode Franck-Condon weighted density of states from the computed S0 Hessian.

The paper's rates use a single-mode MLJ line shape with assumed parameters (lambda_s = 0.10 eV, hbar*omega = 0.17 eV,
S = 0.8). Here the line shape of each S1 -> T_n channel is built from the B3LYP/def2-SVP S0 normal modes
(paper/references/run_s0_hessian.py) in the displaced-oscillator model (equal frequencies in both states):
  * displacements between the two states, projected on the mass-weighted normal modes, either
      (minima)   from the optimized DeltaSCF S1 and T1 minima of run_adiabatic.py, or
      (gradient) from TDA gradients at the S0 minimum (run_vertical_gradients.py), dQ_X = -L^T M^-1/2 g_X / omega^2;
  * Huang-Rhys factors S_k = omega_k dQ_k^2 / 2 (atomic units), reorganization lambda_k = S_k hbar omega_k;
  * modes below 200 cm^-1 and an optional solvent term are treated classically (Marcus Gaussian), the rest
    quantum mechanically at 298 K (Kubo generating function); a minimum Gaussian width of 0.01 eV is kept;
  * rho(dE) = (1/pi) Re int_0^inf dt exp(-i dE t) F(t), dE = E_T - E_S1 (negative downhill, as in psscreen.rates),
    F(t) = exp(sum_k S_k[(n_k+1)(e^{-i w_k t} - 1) + n_k (e^{i w_k t} - 1)] - i lam_cl t - lam_cl kT t^2).
The rate of a channel is k = (2 pi/hbar) V^2 rho(dE) with the pipeline coupling V and gap dE of that channel, and it
is compared with the MLJ rate at the same V and dE.

Usage (repository root): python paper/analysis/fcwd_multimode.py [NAME ...] [--solvent 0.0]
Reads paper/references/results/NAME.s0_hessian.npz, NAME.adiabatic_{S0,S1,T1}.xyz, NAME.vertical_gradients.npz
(whichever exist). Writes paper/analysis/fcwd_multimode.json and .md.
"""
import argparse
import json
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
import numpy as np  # noqa: E402

from psscreen.geometry import read_xyz  # noqa: E402
from psscreen.rates import CM1_TO_EV, HBAR_EVS, KB_EV, MLJParams, fcwd_mlj  # noqa: E402

HA_EV = 27.211386245988
BOHR_A = 0.52917721092
AMU_ME = 1822.888486209
CM1_HA = 1 / 219474.6313705
T = 298.15
KT = KB_EV * T
W_CUT_CM = 200.0
SIGMA_MIN = 0.01
RES = os.path.join(ROOT, "paper", "references", "results")
MASS = {"H": 1.00782503, "B": 11.0093054, "C": 12.0, "N": 14.0030740, "O": 15.9949146, "F": 18.9984032,
        "Cl": 34.9688527, "Br": 78.9183376, "I": 126.904473, "S": 31.9720707}


def normal_modes(H, xyz_A, symbols):
    """Mass-weighted normal modes with translations and rotations projected out. Returns omega (Ha), L (3N x nvib)
    columns in mass-weighted coordinates, sqrt-mass vector."""
    m = np.repeat([MASS[s] * AMU_ME for s in symbols], 3)
    sm = np.sqrt(m)
    Hm = H / np.outer(sm, sm)
    x = np.asarray(xyz_A, float) / BOHR_A
    com = (x * m[0::3, None]).sum(0) / m[0::3].sum()
    r = x - com
    vecs = []
    for k in range(3):
        v = np.zeros_like(x)
        v[:, k] = 1
        vecs.append((v.ravel() * sm))
    for k in range(3):
        e = np.zeros(3)
        e[k] = 1
        vecs.append((np.cross(e, r).ravel() * sm))
    Bq, _ = np.linalg.qr(np.array(vecs).T)
    P = np.eye(len(m)) - Bq @ Bq.T
    w2, L = np.linalg.eigh(P @ Hm @ P)
    keep = np.argsort(np.abs(w2))[6:]
    keep = keep[np.argsort(w2[keep])]
    w2, L = w2[keep], L[:, keep]
    if (w2 < 0).any():
        print(f"  warning: {int((w2 < 0).sum())} imaginary modes, smallest {np.sqrt(-w2.min()) / CM1_HA:.0f}i cm-1;"
              " taken as |omega|", flush=True)
    return np.sqrt(np.abs(w2)), L, sm, P


def kabsch(ref, mob, w):
    """Rotate and translate mob onto ref (mass weights w)."""
    cr, cm = (ref * w[:, None]).sum(0) / w.sum(), (mob * w[:, None]).sum(0) / w.sum()
    A, B = ref - cr, mob - cm
    U, _, Vt = np.linalg.svd((B * w[:, None]).T @ A)
    d = np.sign(np.linalg.det(U @ Vt))
    R = U @ np.diag([1, 1, d]) @ Vt
    return B @ R + cr


def huang_rhys(omega, dQ):
    S = 0.5 * omega * dQ ** 2
    return S, S * omega * HA_EV


def rho_multimode(dE_eV, omega_Ha, S, lam_solv=0.0):
    """FCWD (1/eV) at dE = E_final - E_initial (eV) for displaced oscillators."""
    w = omega_Ha * HA_EV
    lo = w < W_CUT_CM * CM1_TO_EV
    lam_cl = float((S[lo] * w[lo]).sum()) + lam_solv
    lam_cl = max(lam_cl, SIGMA_MIN ** 2 / (2 * KT))          # Gaussian variance 2 lam kT >= sigma_min^2
    wq, Sq = w[~lo], S[~lo]
    n = 1 / np.expm1(wq / KT)
    tmax = 8 / math.sqrt(2 * lam_cl * KT)
    t = np.linspace(0, tmax, 200001)                          # t in 1/eV (hbar = 1)
    ph = np.zeros_like(t, dtype=complex)
    for wk, sk, nk in zip(wq, Sq, n):
        ph += sk * ((nk + 1) * np.expm1(-1j * wk * t) + nk * np.expm1(1j * wk * t))
    F = np.exp(ph - 1j * lam_cl * t - lam_cl * KT * t ** 2)
    dEs = np.atleast_1d(dE_eV)
    out = np.array([np.trapezoid(np.exp(1j * (-d) * t) * F, t).real / math.pi for d in dEs])
    return out if np.ndim(dE_eV) else float(out[0]), lam_cl


def pipeline_channels(name):
    for p in (os.path.join(ROOT, "validation", "runs", f"{name}.json"),
              os.path.join(ROOT, "screening", "runs", "candidates_round1", f"{name}.json")):
        if os.path.exists(p):
            return json.load(open(p))["photophysics"]["channels"]
    raise FileNotFoundError(name)


def analyse(name, lam_solv):
    hz = np.load(os.path.join(RES, f"{name}.s0_hessian.npz"))
    sym = [str(s) for s in hz["symbols"]]
    x0 = np.asarray(hz["xyz_A"], float)
    omega, L, sm, P = normal_modes(hz["H"], x0, sym)
    m_at = np.array([MASS[s] for s in sym])
    res = {"name": name, "n_modes": int(len(omega)), "omega_min_cm": float(omega.min() / CM1_HA),
           "lambda_solvent_eV": lam_solv, "displacements": {}, "channels": []}
    dQ = {}
    # (minima) optimized DeltaSCF minima
    for st in ("S1", "T1"):
        f = os.path.join(RES, f"{name}.adiabatic_{st}.xyz")
        if os.path.exists(f):
            xs = kabsch(x0, np.asarray(read_xyz(f)[1], float), m_at)
            d = (xs - x0).ravel() / BOHR_A * sm
            dQ[("min", st)] = L.T @ (P @ d)
    # (gradient) TDA gradients at the S0 minimum
    fg = os.path.join(RES, f"{name}.vertical_gradients.npz")
    if os.path.exists(fg):
        g = np.load(fg)
        for st in [k for k in g.files if k[0] in "ST" and k[1:].isdigit()]:
            gm = P @ (np.asarray(g[st]).ravel() / sm)
            dQ[("grad", st)] = -(L.T @ gm) / omega ** 2
    for (kind, st), q in dQ.items():
        S, lam = huang_rhys(omega, q)
        res["displacements"][f"{kind}:{st}"] = {"lambda_from_S0_eV": float(lam.sum()),
                                               "S_total": float(S.sum())}
    chans = {c["T"]: c for c in pipeline_channels(name)}
    pairs = [("min", "S1", "T1", 1)] + [("grad", "S1", f"T{n}", n) for n in (1, 2, 3)]
    for kind, a, b, n in pairs:
        if (kind, a) not in dQ or (kind, b) not in dQ or n not in chans:
            continue
        S, lam = huang_rhys(omega, dQ[(kind, a)] - dQ[(kind, b)])
        w_cm = omega / CM1_HA
        top = np.argsort(lam)[::-1][:5]
        c = chans[n]
        dEs = np.linspace(-1.0, 0.5, 301)
        rho_curve, lam_cl = rho_multimode(dEs, omega, S, lam_solv)
        r_mm = float(np.interp(c["dE_eV"], dEs, rho_curve))
        r_mlj = fcwd_mlj(c["dE_eV"], MLJParams())
        V = c["V_cm"] * CM1_TO_EV
        k_mm = 2 * math.pi / HBAR_EVS * V ** 2 * r_mm
        res["channels"].append({
            "model": kind, "channel": f"{a}->{b}", "dE_pipeline_eV": c["dE_eV"], "V_cm": c["V_cm"],
            "lambda_total_eV": float(lam.sum()), "lambda_classical_eV": lam_cl,
            "lambda_quantum_eV": float(lam[w_cm >= W_CUT_CM].sum()),
            "S_quantum": float(S[w_cm >= W_CUT_CM].sum()),
            "mean_quantum_mode_cm": float((lam[w_cm >= W_CUT_CM] @ w_cm[w_cm >= W_CUT_CM]) / max(lam[w_cm >= W_CUT_CM].sum(), 1e-30)),
            "top_modes": [{"omega_cm": float(w_cm[i]), "S": float(S[i]), "lambda_eV": float(lam[i])} for i in top],
            "rho_multimode_per_eV": r_mm, "rho_mlj_per_eV": r_mlj, "k_multimode": k_mm, "k_pipeline_channel": c["k"],
            "ratio_multimode_to_mlj": r_mm / r_mlj if r_mlj > 0 else None,
            "rho_curve": {"dE_eV": dEs[::10].tolist(), "rho": rho_curve[::10].tolist()}})
    return res


def selftest():
    """Single quantum mode at T -> 0 against the analytic MLJ form (classical width lam_s)."""
    p = MLJParams()
    w = np.array([p.hw / HA_EV])
    S = np.array([p.S])
    lam_s = p.lambda_s
    r, _ = rho_multimode(np.array([-0.4, -0.2, 0.0]), w, S, lam_solv=lam_s)
    ref = [fcwd_mlj(d, p) for d in (-0.4, -0.2, 0.0)]
    # MLJ uses n = 0 initial vibrational state; at 298 K the 0.17 eV mode has n = 1e-3, so agreement to ~1e-3
    return [float(x) for x in r], ref


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("names", nargs="*", default=["BDP", "I2-BDP", "azaBDP_Br2"])
    ap.add_argument("--solvent", type=float, default=0.0, help="classical outer-sphere lambda (eV)")
    a = ap.parse_args()
    r, ref = selftest()
    print("self-test (multimode vs MLJ):", [f"{x:.4g}/{y:.4g}" for x, y in zip(r, ref)], flush=True)
    out = {"selftest": {"multimode": r, "mlj": ref}, "molecules": []}
    for n in a.names:
        if not os.path.exists(os.path.join(RES, f"{n}.s0_hessian.npz")):
            print(n, "no Hessian yet", flush=True)
            continue
        out["molecules"].append(analyse(n, a.solvent))
    json.dump(out, open(os.path.join(ROOT, "paper", "analysis", "fcwd_multimode.json"), "w"), indent=1)
    L = ["# Multimode Franck-Condon weighted density of states (generated by fcwd_multimode.py)", "",
         "| Dye | model | channel | dE (eV) | V (cm⁻¹) | λ total (eV) | λ quantum | S quantum | mean ħω (cm⁻¹) | "
         "ρ multimode / MLJ | k multimode (s⁻¹) | k pipeline (s⁻¹) |", "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for m in out["molecules"]:
        for c in m["channels"]:
            L.append(f"| {m['name']} | {c['model']} | {c['channel']} | {c['dE_pipeline_eV']:.3f} | {c['V_cm']:.2f} | "
                     f"{c['lambda_total_eV']:.3f} | {c['lambda_quantum_eV']:.3f} | {c['S_quantum']:.2f} | "
                     f"{c['mean_quantum_mode_cm']:.0f} | {c['ratio_multimode_to_mlj']:.2g} | {c['k_multimode']:.2e} | "
                     f"{c['k_pipeline_channel']:.2e} |")
    open(os.path.join(ROOT, "paper", "analysis", "fcwd_multimode.md"), "w").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
