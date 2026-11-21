"""Sensitivity of k_ISC and Phi_T to the treatment of the upper triplets.

The pipeline replaces E(S1) and E(T1) with Delta-SCF values but keeps the TDDFT (TDA) energies of
T2, T3, ... This script recomputes k_ISC and Phi_T for every dye with two consistent alternatives:

  (a) "T1 shift": every T_n is shifted by the Delta-SCF correction of T1 (E_T1^DSCF - E_T1^TDA),
      i.e. the upper triplets are assumed to carry the same error as T1;
  (b) "S1 shift": every T_n with n >= 2 is shifted by the Delta-SCF correction of S1, i.e. the upper
      triplets are assumed to carry the (large) cyanine-type error of S1.

SOC, k_r and k_ic are kept fixed. Usage (repository root): python paper/analysis/tn_shift_check.py
Outputs: paper/analysis/tn_shift_check.json and .md
"""
from __future__ import annotations
import json, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from psscreen.rates import MLJParams, k_isc  # noqa: E402

RUNS = ROOT / "validation" / "runs"
CAND = ROOT / "screening" / "runs" / "candidates_round1"
DYES = ["BODIPY_core", "BDP", "I-BDP", "I2-BDP", "I3-BDP", "azaBDP_Ph4", "BDP-0", "I-BDP-0", "I2-BDP-0"]
CANDIDATES = ["Cl2-BDP", "Br2-BDP", "I1-BDP", "azaBDP_Br2", "azaBDP_BrPh2", "azaBDP_IPh2"]
P = MLJParams()


def load(name):
    for p in (RUNS / f"{name}.json", RUNS / "set3" / f"{name}_c0.json", CAND / f"{name}.json"):
        if p.exists():
            return json.loads(p.read_text())
    raise FileNotFoundError(name)


def rates(d, shift_upper):
    ph = d["photophysics"]
    k, chans = 0.0, []
    for c in ph["channels"]:
        dE = c["dE_eV"] + (shift_upper if c["T"] >= 2 else 0.0)
        kc = k_isc(c["V_cm"], dE, P)
        k += kc
        chans.append((c["T"], dE, kc))
    phi = k / (k + ph["k_r"] + ph["k_ic"])
    dom = max(chans, key=lambda x: x[2])
    return k, phi, dom


def main():
    rows = []
    for name in DYES + CANDIDATES:
        d = load(name)
        corr = d["corrections"]
        cS, cT = corr["S"], corr["T"]
        # in the pipeline, T1 carries cT and T_n>=2 carry no correction; S1 carries cS
        # (a): T_n>=2 also get cT  ->  shift of upper channels = +cT
        # (b): T_n>=2 get cS       ->  shift of upper channels = +cS
        base = rates(d, 0.0)
        a = rates(d, cT)
        b = rates(d, cS)
        rows.append({"dye": name, "corr_S1_eV": cS, "corr_T1_eV": cT,
                     "pipeline": {"k": base[0], "phi": base[1], "dom_T": base[2][0], "dom_dE": base[2][1]},
                     "T1_shift": {"k": a[0], "phi": a[1], "dom_T": a[2][0], "dom_dE": a[2][1]},
                     "S1_shift": {"k": b[0], "phi": b[1], "dom_T": b[2][0], "dom_dE": b[2][1]}})
    out = ROOT / "paper" / "analysis"
    (out / "tn_shift_check.json").write_text(json.dumps(rows, indent=1))
    L = ["# Upper-triplet treatment: sensitivity of k_ISC and Phi_T", "",
         "Pipeline: T_n>=2 keep TDDFT energies. (a) T_n>=2 shifted by the Delta-SCF correction of T1. "
         "(b) T_n>=2 shifted by the Delta-SCF correction of S1. Rates in s^-1; dominant channel T_n and its gap (eV).", "",
         "| Dye | corr S1 | corr T1 | k pipeline | Phi_T | dom. | k (a) | Phi_T (a) | dom. (a) | k (b) | Phi_T (b) | dom. (b) |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        f = lambda s: f"{r[s]['k']:.2e} | {r[s]['phi']:.2f} | T{r[s]['dom_T']} ({r[s]['dom_dE']:+.2f})"
        L.append(f"| {r['dye']} | {r['corr_S1_eV']:+.3f} | {r['corr_T1_eV']:+.3f} | {f('pipeline')} | {f('T1_shift')} | {f('S1_shift')} |")
    (out / "tn_shift_check.md").write_text("\n".join(L) + "\n")
    lab = {"BODIPY_core": "BODIPY core", "azaBDP_Ph4": "aza-BDP-Ph$_4$", "azaBDP_Br2": "aza-BDP-Br$_2$",
           "azaBDP_BrPh2": "aza-BDP-BrPh$_2$", "azaBDP_IPh2": "aza-BDP-IPh$_2$"}

    def sci(x):
        m, e = f"{x:.1e}".split("e")
        return f"${m}\\times10^{{{int(e)}}}$"
    T = [r"\section{Treatment of the upper triplets}",
         r"\begin{table}",
         r"  \caption{$\kisc$ (s$^{-1}$) and $\Phi_T$ of the pipeline (T$_{n\ge2}$ at TDA energies) and with all T$_{n\ge2}$",
         r"  shifted by the $\Delta$SCF correction of T$_1$ (a) or of S$_1$ (b); corrections in eV. Dominant channel and its",
         r"  gap (eV) in parentheses. Script \texttt{paper/analysis/tn\_shift\_check.py}.}",
         r"  \label{tab:tnshift}",
         r"  \scriptsize\setlength{\tabcolsep}{2pt}",
         r"  \begin{tabular}{lrrlllll}",
         r"    \toprule",
         r"    Dye & corr.\ S$_1$ & corr.\ T$_1$ & pipeline & (a) T$_1$ shift & (b) S$_1$ shift & $\Phi_T$ pipe./(a)/(b) \\",
         r"    \midrule"]
    for r in rows:
        c = lambda s: f"{sci(r[s]['k'])} (T$_{r[s]['dom_T']}$, ${r[s]['dom_dE']:+.2f}$)"
        T.append(f"    {lab.get(r['dye'], r['dye'])} & ${r['corr_S1_eV']:+.2f}$ & ${r['corr_T1_eV']:+.2f}$ & {c('pipeline')} & "
                 f"{c('T1_shift')} & {c('S1_shift')} & {r['pipeline']['phi']:.2f}/{r['T1_shift']['phi']:.2f}/{r['S1_shift']['phi']:.2f} \\\\")
    T += [r"    \bottomrule", r"  \end{tabular}", r"\end{table}", ""]
    (ROOT / "paper" / "si_tn_shift.tex").write_text("\n".join(T) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
