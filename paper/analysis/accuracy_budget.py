"""Accuracy budget: how accurate must excited-state energies be for photosensitizer screening?

For every validation dye this script takes the computed S1->Tn channels (gap dE_n and SOC V_n) and the
radiative/IC rates, perturbs the energies, and recomputes k_ISC and Phi_T with the same
Marcus-Levich-Jortner (MLJ) model as the pipeline. It reports

  A. the factor by which k_ISC changes for a +/-0.1 eV error in the S1-T gaps, and the largest gap error
     that keeps k_ISC within a factor of 2 and of 10 (two error models: all triplets shifted together
     relative to S1, or only the upper triplets T2, T3, ... shifted);
  B. the largest gap error that keeps the triplet-yield class (Phi_T above or below 0.5) unchanged;
  C. the distance of S1 from the edges of the 650-850 nm therapeutic window (the absorption margin);
  D. the phase-estimation precision per state that these tolerances imply, and the Toffoli cost that
     follows from scaling our CAS(12,11) double-factorization estimate (4.4e7 Toffolis at 1.6 mHa) as 1/eps.

Usage (from the repository root):
    python paper/analysis/accuracy_budget.py [--xc-dir DIR]
DIR is a folder with per-functional subfolders (pbe0/, bhandhlyp/) as in validation/runs_xc on the
results-tars branch; if given, the absorption margin and k_ISC are also reported for those functionals.
Outputs: paper/analysis/accuracy_budget.json, paper/analysis/accuracy_budget.md and paper/figures/*.pdf.
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from psscreen.rates import MLJParams, k_isc  # noqa: E402

RUNS = ROOT / "validation" / "runs"
OUT = ROOT / "paper" / "analysis"
FIG = ROOT / "paper" / "figures"
EV_NM = 1239.84198
HA_EV = 27.211386245988
REF = json.loads((ROOT / "validation" / "references.json").read_text())

DYES = ["BODIPY_core", "BDP", "I-BDP", "I2-BDP", "I3-BDP", "azaBDP_Ph4", "BDP-0", "I-BDP-0", "I2-BDP-0"]
CAND_DIR = ROOT / "screening" / "runs" / "candidates_round1"
CANDIDATES = ["Cl2-BDP", "Br2-BDP", "I1-BDP", "azaBDP_Br2", "azaBDP_BrPh2", "azaBDP_IPh2"]
WINDOW_EV = (EV_NM / 850.0, EV_NM / 650.0)          # 1.459 - 1.907 eV
QPE_REF = {"eps_Ha": 1.6e-3, "toffoli": 4.4e7}     # CAS(12,11), I2-BDP, docs/VALIDATION.md
MLJ_SETS = {"lambda_s=0.05": MLJParams(lambda_s=0.05), "lambda_s=0.10 (default)": MLJParams(),
            "lambda_s=0.20": MLJParams(lambda_s=0.20)}


def load(name, base=RUNS):
    for p in (base / f"{name}.json", base / "set3" / f"{name}_c0.json"):
        if p.exists():
            return json.loads(p.read_text())
    return None


def rates(d, delta_all=0.0, delta_upper=0.0, p=MLJParams()):
    """k_ISC and Phi_T after shifting all triplets (delta_all) and/or T2 and above (delta_upper), in eV,
    relative to S1. Positive delta moves triplets up."""
    ph = d["photophysics"]
    k = 0.0
    for c in ph["channels"]:
        dE = c["dE_eV"] + delta_all + (delta_upper if c["T"] >= 2 else 0.0)
        k += k_isc(c["V_cm"], dE, p)
    phi = k / (k + ph["k_r"] + ph["k_ic"])
    return k, phi


def tolerance(d, factor, mode, p, dmax=0.5, n=501):
    """Largest |delta| (eV) such that k_ISC stays within `factor` of its unperturbed value for every
    error of that size or smaller, in either direction."""
    k0, _ = rates(d, p=p)
    grid = np.linspace(0, dmax, n)
    for sign in (+1, -1):
        for x in grid:
            kw = {"delta_all": sign * x} if mode == "all" else {"delta_upper": sign * x}
            k, _ = rates(d, p=p, **kw)
            if k0 > 0 and (k / k0 > factor or k0 / max(k, 1e-300) > factor):
                grid = grid[grid < x]
                break
    return float(grid[-1]) if len(grid) else 0.0


def phi_margin(d, mode, p, dmax=0.5, n=501):
    """Largest |delta| that keeps Phi_T on the same side of 0.5."""
    _, phi0 = rates(d, p=p)
    side = phi0 > 0.5
    best = dmax
    for sign in (+1, -1):
        for x in np.linspace(0, dmax, n):
            kw = {"delta_all": sign * x} if mode == "all" else {"delta_upper": sign * x}
            _, phi = rates(d, p=p, **kw)
            if (phi > 0.5) != side:
                best = min(best, x)
                break
    return float(best)


def window_margin(E):
    lo, hi = WINDOW_EV
    inside = lo <= E <= hi
    return inside, float(min(abs(E - lo), abs(E - hi)))


def qpe_cost(gap_tol_eV):
    """Per-state QPE precision if the two state energies have independent errors adding in quadrature,
    and the Toffoli count from the 1/eps scaling of the reference estimate."""
    eps_Ha = gap_tol_eV / np.sqrt(2) / HA_EV
    return eps_Ha, QPE_REF["toffoli"] * QPE_REF["eps_Ha"] / eps_Ha


def analyse(base=RUNS, names=DYES):
    rows = []
    for n in names:
        d = load(n, base)
        if d is None:
            continue
        k0, phi0 = rates(d)
        r = {"name": n, "E_S1_eV": d["E_S1_eV"],
             "E_exp_eV": EV_NM / REF[n]["lambda_abs_nm"] if n in REF else None,
             "k_isc": k0, "phi_T": phi0,
             "dominant": max(d["photophysics"]["channels"], key=lambda c: c["k"])}
        for mode in ("all", "upper"):
            kp, _ = rates(d, **({"delta_all": 0.1} if mode == "all" else {"delta_upper": 0.1}))
            km, _ = rates(d, **({"delta_all": -0.1} if mode == "all" else {"delta_upper": -0.1}))
            r[f"factor_+0.1_{mode}"] = kp / k0
            r[f"factor_-0.1_{mode}"] = km / k0
            r[f"phi_+0.1_{mode}"] = rates(d, **({"delta_all": 0.1} if mode == "all" else {"delta_upper": 0.1}))[1]
            r[f"phi_-0.1_{mode}"] = rates(d, **({"delta_all": -0.1} if mode == "all" else {"delta_upper": -0.1}))[1]
            for lab, p in MLJ_SETS.items():
                r.setdefault("tol2", {}).setdefault(mode, {})[lab] = tolerance(d, 2.0, mode, p)
                r.setdefault("tol10", {}).setdefault(mode, {})[lab] = tolerance(d, 10.0, mode, p)
                r.setdefault("phi_margin", {}).setdefault(mode, {})[lab] = phi_margin(d, mode, p)
        inside, m = window_margin(d["E_S1_eV"])
        r["window_inside"], r["window_margin_eV"] = inside, m
        if r["E_exp_eV"] is not None:
            insx, mx = window_margin(r["E_exp_eV"])
            r["window_inside_exp"], r["window_margin_exp_eV"] = insx, mx
        rows.append(r)
    return rows


def fig_fcwd(rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dE = np.linspace(-1.2, 0.8, 801)
    fig, ax = plt.subplots(figsize=(3.4, 3.0))
    for (lab, p), col in zip(MLJ_SETS.items(), ["#9f7aea", "#2b6cb0", "#38a169"]):
        y = [k_isc(1.0, x, p) for x in dE]
        ax.semilogy(dE, y, color=col, lw=1.2, label=lab.replace("lambda_s", r"$\lambda_s$"))
    marks = {"I2-BDP": ("o", "#c05621"), "I2-BDP-0": ("s", "#c05621")}
    for r in rows:
        if r["name"] not in marks:
            continue
        mk, col = marks[r["name"]]
        d = load(r["name"])
        for c in d["photophysics"]["channels"][:3]:
            if abs(c["dE_eV"]) < 1.2:
                ax.plot(c["dE_eV"], k_isc(1.0, c["dE_eV"]), mk, color=col, ms=4,
                        mfc="none" if r["name"].endswith("-0") else col)
                ax.annotate(f"T{c['T']}", (c["dE_eV"], k_isc(1.0, c["dE_eV"])), fontsize=6,
                            textcoords="offset points",
                            xytext=(-12, -9) if r["name"].endswith("-0") else (3, 3))
    ax.set_ylim(1e-12, 1e10)
    ax.set_xlabel(r"$E(\mathrm{T}_n)-E(\mathrm{S}_1)$ (eV)")
    ax.set_ylabel(r"$k_\mathrm{ISC}$ for $V$ = 1 cm$^{-1}$ (s$^{-1}$)")
    ax.legend(frameon=False, fontsize=6.5, loc="lower left")
    ax.text(0.98, 0.97, "filled: I2-BDP\nopen: I2-BDP-0", transform=ax.transAxes, ha="right", va="top",
            fontsize=6)
    fig.tight_layout()
    fig.savefig(FIG / "mlj_gap_sensitivity.pdf")


def fmt_m(x, dmax=0.5):
    return f">{dmax:.1f}" if x >= dmax else f"{x:.2f}"


def markdown(rows, xc_rows):
    L = ["# Accuracy budget (generated by paper/analysis/accuracy_budget.py)", "",
         "Error model `all`: every triplet shifted by delta relative to S1 (S-T gap error). "
         "Model `upper`: only T2, T3, ... shifted (error in the upper triplets). Default MLJ parameters "
         "unless stated; tolerances in eV.", "",
         "## A/B. k_ISC and Phi_T sensitivity (B3LYP, default MLJ)", "",
         "| Dye | k_ISC (s⁻¹) | Φ_T | dominant channel | k(+0.1)/k, k(−0.1)/k [all] | tol ×2 [all] | tol ×10 [all] | "
         "k(+0.1)/k, k(−0.1)/k [upper] | tol ×2 [upper] | Φ_T class margin [all / upper] |",
         "|---|---|---|---|---|---|---|---|---|---|"]
    D = "lambda_s=0.10 (default)"
    for r in rows:
        if r["name"] == CANDIDATES[0]:
            L.append("| *screening candidates (no experiment)* | | | | | | | | | |")
        dm = r["dominant"]
        L.append(f"| {r['name']} | {r['k_isc']:.2g} | {r['phi_T']:.2f} | T{dm['T']} ({dm['dE_eV']:+.2f} eV) | "
                 f"{r['factor_+0.1_all']:.2g}, {r['factor_-0.1_all']:.2g} | {r['tol2']['all'][D]:.3f} | "
                 f"{r['tol10']['all'][D]:.3f} | {r['factor_+0.1_upper']:.2g}, {r['factor_-0.1_upper']:.2g} | "
                 f"{r['tol2']['upper'][D]:.3f} | {fmt_m(r['phi_margin']['all'][D])} / {fmt_m(r['phi_margin']['upper'][D])} |")
    L += ["", "## Dependence of the ×2 tolerance on λ_s (error model `all`)", "",
          "| Dye | " + " | ".join(MLJ_SETS) + " |", "|---|" + "---|" * len(MLJ_SETS)]
    for r in rows:
        L.append(f"| {r['name']} | " + " | ".join(f"{r['tol2']['all'][k]:.3f}" for k in MLJ_SETS) + " |")
    L += ["", "## C. Absorption-window margin (650-850 nm = 1.459-1.907 eV)", "",
          "| Dye | E(S1) calc | inside (calc) | margin calc (eV) | E exp | inside (exp) | margin exp (eV) |",
          "|---|---|---|---|---|---|---|"]
    for r in rows:
        if r["E_exp_eV"] is None:
            continue
        L.append(f"| {r['name']} | {r['E_S1_eV']:.2f} | {r['window_inside']} | {r['window_margin_eV']:.2f} | "
                 f"{r['E_exp_eV']:.2f} | {r['window_inside_exp']} | {r['window_margin_exp_eV']:.2f} |")
    # D: QPE
    heavy = [r for r in rows if r["k_isc"] > 1e6]
    L += ["", "## D. Implied phase-estimation precision (from the ×2 k_ISC tolerance, model `all`)", "",
          "Per-state error eps = tol/√2 (independent errors on S1 and T_n); Toffolis scaled as 1/eps from "
          f"{QPE_REF['toffoli']:.1e} at {QPE_REF['eps_Ha']*1e3:.1f} mHa (CAS(12,11), I2-BDP).", "",
          "| Dye | gap tol ×2 (eV) | eps per state (mHa) | Toffolis per QPE run |", "|---|---|---|---|"]
    for r in heavy:
        e, t = qpe_cost(r["tol2"]["all"][D])
        L.append(f"| {r['name']} | {r['tol2']['all'][D]:.3f} | {e*1e3:.2f} | {t:.1e} |")
    if xc_rows:
        L += ["", "## Functional dependence (ΔSCF-corrected S1 error and k_ISC)", "",
              "| Dye | " + " | ".join(f"{x} err (eV) | {x} k_ISC" for x in xc_rows) + " | k_ISC exp |",
              "|---|" + "---|---|" * len(xc_rows) + "---|"]
        for n in DYES:
            cells = []
            for x, rr in xc_rows.items():
                m = {r["name"]: r for r in rr}.get(n)
                cells += [f"{m['E_S1_eV']-m['E_exp_eV']:+.2f}", f"{m['k_isc']:.1e}"] if m else ["—", "—"]
            kexp = REF[n].get("k_isc_s")
            kexp = f"{kexp:.1e}" if isinstance(kexp, (int, float)) else (kexp or "—")
            L.append(f"| {n} | " + " | ".join(cells) + f" | {kexp} |")
    return "\n".join(L) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xc-dir", type=Path, default=None)
    a = ap.parse_args()
    rows = analyse() + analyse(CAND_DIR, CANDIDATES)
    xc_rows = {"B3LYP": [r for r in rows if r["name"] in DYES]}
    if a.xc_dir:
        for x in ("pbe0", "bhandhlyp"):
            if (a.xc_dir / x).exists():
                xc_rows[x.upper()] = analyse(a.xc_dir / x)
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    (OUT / "accuracy_budget.json").write_text(json.dumps({"B3LYP": rows, **{k: v for k, v in xc_rows.items()
                                                                          if k != "B3LYP"}}, indent=1, default=float))
    (OUT / "accuracy_budget.md").write_text(markdown(rows, xc_rows if len(xc_rows) > 1 else None))
    fig_fcwd(rows)
    print((OUT / "accuracy_budget.md").read_text())


if __name__ == "__main__":
    main()
