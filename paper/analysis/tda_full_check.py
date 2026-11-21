"""Effect of the truncated TDA window on the upper triplets, k_ISC and Phi_T.

The pipeline obtains T_n (n >= 2) from a windowed TDA (10 eV window, at most 40 occupied x 60 virtual orbitals).
paper/references/run_tda_full.py recomputes the triplets of selected dyes with the same window (control), a 20 eV
window without caps, and the full TDA. This script shifts every channel T_n >= 2 of the pipeline by the
full-minus-window difference of that state and recomputes k_ISC and Phi_T (SOC, k_r, k_ic fixed).
Usage (repository root): python paper/analysis/tda_full_check.py
Outputs: paper/analysis/tda_full_check.json and .md, paper/si_tda_full.tex
"""
from __future__ import annotations
import json, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from psscreen.rates import MLJParams, k_isc  # noqa: E402

RES = ROOT / "paper" / "references" / "results"
RUNS = ROOT / "validation" / "runs"
CAND = ROOT / "screening" / "runs" / "candidates_round1"
P = MLJParams()
ORDER = ["BODIPY_core", "I2-BDP", "I2-BDP-0", "azaBDP_Br2"]
LAB = {"BODIPY_core": "BODIPY core", "azaBDP_Br2": "aza-BDP-Br$_2$"}


def load_run(name):
    for p in (RUNS / f"{name}.json", RUNS / "set3" / f"{name}_c0.json", CAND / f"{name}.json"):
        if p.exists():
            return json.loads(p.read_text())
    raise FileNotFoundError(name)


def rates(d, shifts):
    ph = d["photophysics"]
    k, best = 0.0, None
    for c in ph["channels"]:
        dE = c["dE_eV"] + (shifts.get(c["T"], 0.0) if c["T"] >= 2 else 0.0)
        kc = k_isc(c["V_cm"], dE, P)
        k += kc
        if best is None or kc > best[2]:
            best = (c["T"], dE, kc)
    return k, k / (k + ph["k_r"] + ph["k_ic"]), best


def main():
    rows = []
    for name in ORDER:
        f = RES / f"{name}.tda_full.json"
        if not f.exists():
            continue
        t = json.loads(f.read_text())
        win, w20, full = t["window10"]["E_T_eV"], t["window20"]["E_T_eV"], t["full"]["E_T_eV"]
        n = min(len(win), len(full), 5)
        dfull = {i + 1: full[i] - win[i] for i in range(n)}
        d20 = {i + 1: w20[i] - win[i] for i in range(n)}
        ctrl = max(abs(a - b) for a, b in zip(t["pipeline_E_T_eV"], win))
        d = load_run(name)
        base, new = rates(d, {}), rates(d, dfull)
        chans = {c["T"]: c["dE_eV"] for c in d["photophysics"]["channels"]}
        rows.append({"dye": name, "control_max_abs_eV": ctrl, "window": t.get("pipeline_window"),
                     "dE_20_minus_10_eV": d20, "dE_full_minus_10_eV": dfull,
                     "gap_T2_S1_pipeline_eV": chans.get(2), "gap_T2_S1_full_eV": chans.get(2, 0) + dfull[2],
                     "pipeline": {"k": base[0], "phi": base[1], "dom_T": base[2][0], "dom_dE": base[2][1]},
                     "full": {"k": new[0], "phi": new[1], "dom_T": new[2][0], "dom_dE": new[2][1]}})
    out = ROOT / "paper" / "analysis"
    (out / "tda_full_check.json").write_text(json.dumps(rows, indent=1))
    L = ["# Truncated TDA window: effect on T_n, k_ISC and Phi_T", "",
         "| Dye | dT2 (20 eV / full) | dT3 (20 eV / full) | T2-S1 pipeline -> full | k pipeline | k full-TDA | ratio | Phi_T pipe -> full |",
         "|---|---|---|---|---|---|---|---|"]
    for r in rows:
        a, b = r["dE_20_minus_10_eV"], r["dE_full_minus_10_eV"]
        L.append(f"| {r['dye']} | {a[2]:+.3f} / {b[2]:+.3f} | {a[3]:+.3f} / {b[3]:+.3f} | "
                 f"{r['gap_T2_S1_pipeline_eV']:+.3f} -> {r['gap_T2_S1_full_eV']:+.3f} | {r['pipeline']['k']:.2e} | "
                 f"{r['full']['k']:.2e} | {r['full']['k'] / r['pipeline']['k']:.2f} | {r['pipeline']['phi']:.2f} -> {r['full']['phi']:.2f} |")
    (out / "tda_full_check.md").write_text("\n".join(L) + "\n")

    def sci(x):
        m, e = f"{x:.1e}".split("e")
        return f"${m}\\times10^{{{int(e)}}}$"
    T = [r"\section{Convergence of the triplet TDA window}",
         r"\begin{table}",
         r"  \caption{Shift of T$_2$ and T$_3$ (eV) when the pipeline TDA window (10~eV, at most 40 occupied $\times$ 60",
         r"  virtual orbitals) is widened to 20~eV without caps and replaced by the full TDA (B3LYP/def2-SVP, same geometry),",
         r"  the resulting T$_2$--S$_1$ gap (eV), and $\kisc$ (s$^{-1}$) and $\Phi_T$ with every T$_{n\ge2}$ shifted by its",
         r"  full-TDA correction. Scripts \texttt{paper/references/run\_tda\_full.py} and \texttt{paper/analysis/tda\_full\_check.py}.}",
         r"  \label{tab:tdafull}",
         r"  \footnotesize\setlength{\tabcolsep}{3pt}",
         r"  \begin{tabular}{lrrrrrrll}",
         r"    \toprule",
         r"    & \multicolumn{2}{c}{$\Delta$T$_2$} & \multicolumn{2}{c}{$\Delta$T$_3$} & \multicolumn{2}{c}{T$_2-$S$_1$} & \multicolumn{2}{c}{$\kisc$ ($\Phi_T$)} \\",
         r"    Dye & 20~eV & full & 20~eV & full & pipeline & full & pipeline & full TDA \\",
         r"    \midrule"]
    for r in rows:
        a, b = r["dE_20_minus_10_eV"], r["dE_full_minus_10_eV"]
        T.append(f"    {LAB.get(r['dye'], r['dye'])} & ${a[2]:+.3f}$ & ${b[2]:+.3f}$ & ${a[3]:+.3f}$ & ${b[3]:+.3f}$ & "
                 f"${r['gap_T2_S1_pipeline_eV']:+.3f}$ & ${r['gap_T2_S1_full_eV']:+.3f}$ & "
                 f"{sci(r['pipeline']['k'])} ({r['pipeline']['phi']:.2f}) & {sci(r['full']['k'])} ({r['full']['phi']:.2f}) \\\\")
    T += [r"    \bottomrule", r"  \end{tabular}", r"\end{table}", ""]
    (ROOT / "paper" / "si_tda_full.tex").write_text("\n".join(T) + "\n")
    print("\n".join(L))
    for r in rows:
        print(r["dye"], "control max |pipeline - window10| =", f"{r['control_max_abs_eV']:.1e}", r["window"])


if __name__ == "__main__":
    main()
