"""SI table: functional (B3LYP, PBE0, BHandHLYP; Delta-SCF S1/T1 in each) and conformer dependence.

Reads validation/runs (B3LYP), validation/runs_xc/<xc> (from the tars XC scan, validation/run_xc_scan.py) and the
three xTB conformers of the distyryl dyes in validation/runs/set3 (validation/conformers.py).
Usage (repository root): python paper/analysis/xc_conformer_table.py
Output: paper/si_xc_conf.tex
"""
from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
V = ROOT / "validation"
EV = 1239.84198
DYES = ["BODIPY_core", "BDP", "I-BDP", "I2-BDP", "I3-BDP", "azaBDP_Ph4", "BDP-0", "I-BDP-0", "I2-BDP-0"]
LAB = {"BODIPY_core": "BODIPY core", "azaBDP_Ph4": "aza-BDP-Ph$_4$"}
REF = json.loads((V / "references.json").read_text())


def b3lyp(n):
    for p in (V / "runs" / f"{n}.json", V / "runs" / "set3" / f"{n}_c0.json"):
        if p.exists():
            return json.loads(p.read_text())


def sci(x):
    m, e = f"{x:.1e}".split("e")
    return f"${m}\\times10^{{{int(e)}}}$"


def main():
    T = [r"\section{Functional and conformer dependence}",
         r"\begin{table}",
         r"  \caption{Dependence on the exchange--correlation functional (def2-SVP, $\Delta$SCF S$_1$ and T$_1$ with each",
         r"  functional, same xTB geometries): error of S$_1$ against the absorption maximum (eV), S$_1-$T$_1$ gap (eV),",
         r"  $\kisc$ (s$^{-1}$) and $\Phi_T$. Measured $\kisc$: I2-BDP $3.5\times10^{9}$, I3-BDP $3.9\times10^{9}$,",
         r"  I2-BDP-0 $4.5\times10^{8}$~s$^{-1}$. Scripts \texttt{validation/run\_xc\_scan.py} and",
         r"  \texttt{paper/analysis/xc\_conformer\_table.py}.}",
         r"  \label{tab:xc}",
         r"  \scriptsize\setlength{\tabcolsep}{2pt}",
         r"  \begin{tabular}{l" + "rrrr" * 3 + "}",
         r"    \toprule",
         r"    & \multicolumn{4}{c}{B3LYP (20\%)} & \multicolumn{4}{c}{PBE0 (25\%)} & \multicolumn{4}{c}{BHandHLYP (50\%)} \\",
         r"    Dye" + r" & $\delta$S$_1$ & S$_1-$T$_1$ & $\kisc$ & $\Phi_T$" * 3 + r" \\",
         r"    \midrule"]
    for n in DYES:
        cells = []
        for xc in ("b3lyp", "pbe0", "bhandhlyp"):
            r = b3lyp(n) if xc == "b3lyp" else json.loads((V / "runs_xc" / xc / f"{n}.json").read_text())
            ph = r["photophysics"]
            cells.append(f"${r['E_S1_eV'] - EV / REF[n]['lambda_abs_nm']:+.2f}$ & {r['dE_ST_eV']:.2f} & "
                         f"{sci(ph['k_isc'])} & {ph['phi_ISC']:.2f}")
        T.append(f"    {LAB.get(n, n)} & " + " & ".join(cells) + r" \\")
    T += [r"    \bottomrule", r"  \end{tabular}", r"\end{table}", "",
          r"\begin{table}",
          r"  \caption{xTB conformers of the distyryl dyes (B3LYP/$\Delta$SCF): S$_1$ (eV), $\kisc$ (s$^{-1}$) and $\Phi_T$.",
          r"  Conformer c0 is the one used in the main text. Script",
          r"  \texttt{validation/conformers.py}.}",
          r"  \label{tab:conf}",
          r"  \small",
          r"  \begin{tabular}{llrrr}",
          r"    \toprule",
          r"    Dye & conformer & S$_1$ & $\kisc$ & $\Phi_T$ \\",
          r"    \midrule"]
    for n in ("BDP-0", "I-BDP-0", "I2-BDP-0"):
        for c in range(3):
            p = V / "runs" / "set3" / f"{n}_c{c}.json"
            if not p.exists():
                continue
            r = json.loads(p.read_text())
            ph = r["photophysics"]
            T.append(f"    {n} & c{c} & {r['E_S1_eV']:.2f} & {sci(ph['k_isc'])} & {ph['phi_ISC']:.2f} \\\\")
    T += [r"    \bottomrule", r"  \end{tabular}", r"\end{table}", ""]
    (ROOT / "paper" / "si_xc_conf.tex").write_text("\n".join(T) + "\n")
    print("\n".join(T))


if __name__ == "__main__":
    main()
