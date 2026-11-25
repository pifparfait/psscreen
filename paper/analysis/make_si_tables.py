"""Supporting-information tables for all fifteen dyes: pipeline energies, couplings and rates; accuracy budget;
decision costs; and Cartesian coordinates. Writes paper/si_tables.tex (input by si.tex).
Usage (repository root): python paper/analysis/make_si_tables.py
"""
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "paper" / "analysis"))
from accuracy_budget import CAND_DIR, CANDIDATES, DYES, load  # noqa: E402

LABEL = {"BODIPY_core": "BODIPY core", "azaBDP_Ph4": "aza-BDP-Ph$_4$", "azaBDP_Br2": "aza-BDP-Br$_2$",
         "azaBDP_BrPh2": "aza-BDP-BrPh$_2$", "azaBDP_IPh2": "aza-BDP-IPh$_2$"}
ALL = DYES + CANDIDATES


def lab(n):
    return LABEL.get(n, n)


def sci(x, nd=1):
    if x is None:
        return "--"
    if x == 0:
        return "0"
    e = int(math.floor(math.log10(abs(x))))
    m = x / 10 ** e
    if round(m, nd) >= 10:
        m, e = m / 10, e + 1
    return f"${m:.{nd}f}\\times10^{{{e}}}$"


def data(n):
    return load(n) or load(n, CAND_DIR)


def geo_path(n):
    for p in (ROOT / "validation" / "runs" / f"{n}.xyz", ROOT / "validation" / "runs" / "set3" / f"{n}_c0.xyz",
              CAND_DIR / f"{n}.xyz"):
        if p.exists():
            return p
    return None


out = []
w = out.append

# --- pipeline results -----------------------------------------------------------------------------------------
w(r"""\section{Pipeline results for all dyes}
Table~\ref{tab:allpipe} lists the B3LYP/$\Delta$SCF pipeline output for the nine dyes of the benchmark and the six
candidates of the first screening round. The machine-readable files are in the repository, in
\texttt{validation/runs} and \texttt{screening/runs/candidates\_round1}.

\begin{table}
  \caption{B3LYP/def2-SVP pipeline results: $\Delta$SCF S$_1$ and T$_1$, TDA T$_2$ and T$_3$ (eV; raw TDA, before the shift by the $\Delta$SCF correction of T$_1$ that the rates use), S$_1$--T$_n$
  spin--orbit couplings (cm$^{-1}$), ISC rate (s$^{-1}$) and triplet yield. Dyes below the line are the screening
  candidates without experimental data.}
  \label{tab:allpipe}
  \small
  \begin{tabular}{lccccrrrcc}
    \toprule
    Dye & S$_1$ & T$_1$ & T$_2$ & T$_3$ & $V_1$ & $V_2$ & $V_3$ & $\kisc$ & $\Phi_T$ \\
    \midrule""")
for n in ALL:
    if n == CANDIDATES[0]:
        w(r"    \midrule")
    d = data(n)
    T = d["E_triplets_eV"]
    V = d["SOC_S1_Tn_cm"]
    ph = d["photophysics"]
    w(f"    {lab(n)} & {d['E_S1_eV']:.2f} & {d['E_T1_eV']:.2f} & {T[1]:.2f} & {T[2]:.2f} & {V[0]:.2f} & "
      f"{V[1]:.2f} & {V[2]:.2f} & {sci(ph['k_isc'])} & {ph.get('phi_ISC', ph.get('phi_T', 0)):.2f} \\\\")
w(r"""    \bottomrule
  \end{tabular}
\end{table}
""")

# --- accuracy budget ------------------------------------------------------------------------------------------
ab = {r["name"]: r for r in json.loads((ROOT / "paper" / "analysis" / "accuracy_budget.json").read_text())["B3LYP"]}
w(r"""\section{Accuracy budget for all dyes}
\begin{table}
  \caption{Accuracy budget (B3LYP pipeline as reference). Dominant channel T$_n$ and its gap
  $\Delta E=E(\mathrm{T}_n)-E(\mathrm{S}_1)$ (eV); largest rigid shift of all triplets relative to S$_1$ (eV) that keeps
  $\kisc$ within a factor 2 for three values of the reorganization energy $\lambda_s$ (eV), and within a factor 10
  at $\lambda_s=0.10$~eV; factor by which $\kisc$ changes for shifts of $\pm0.1$~eV.}
  \label{tab:allbudget}
  \small
  \begin{tabular}{lccccccc}
    \toprule
    & & & \multicolumn{3}{c}{factor 2} & factor 10 & \\
    \cmidrule(lr){4-6}
    Dye & T$_n$ & $\Delta E$ & $\lambda_s$ 0.05 & 0.10 & 0.20 & 0.10 & $\kisc(\pm0.1)/\kisc$ \\
    \midrule""")
for n in ALL:
    if n == CANDIDATES[0]:
        w(r"    \midrule")
    r = ab[n]
    t2, t10 = r["tol2"]["all"], r["tol10"]["all"]
    w(f"    {lab(n)} & T$_{r['dominant']['T']}$ & ${r['dominant']['dE_eV']:+.2f}$ & {t2['lambda_s=0.05']:.3f} & "
      f"{t2['lambda_s=0.10 (default)']:.3f} & {t2['lambda_s=0.20']:.3f} & {t10['lambda_s=0.10 (default)']:.3f} & "
      f"{r['factor_+0.1_all']:.2g} / {r['factor_-0.1_all']:.2g} \\\\")
w(r"""    \bottomrule
  \end{tabular}
\end{table}
""")

# --- decision costs -------------------------------------------------------------------------------------------
dc = {r["name"]: r for r in json.loads((ROOT / "paper" / "quantum" / "results" / "decision_cost.json").read_text())["rows"]}
w(r"""\section{Decision costs for all dyes}
\begin{table}
  \caption{Decision-level costs (Toffolis) with the model of Section~S4 (THC rank 33, confidence 0.95,
  $\gamma^2=0.5$). $p$: success probability of the filter call; D1: $\kisc$ to 20\%, by amplitude estimation
  (MLAE), by canonical amplitude estimation, and by the channel-resolved route; D2: whether $\Phi_T>0.5$, by MLAE (Table 4 of the main text lists the cheapest estimator, which is sampling for some light-atom dyes).}
  \label{tab:alldecision}
  \footnotesize
  \setlength{\tabcolsep}{4pt}
  \begin{tabular}{lccccccc}
    \toprule
    & & & \multicolumn{3}{c}{D1 ($\kisc$ to 20\%)} & \multicolumn{2}{c}{D2 ($\Phi_T>0.5$)} \\
    \cmidrule(lr){4-6}\cmidrule(lr){7-8}
    Dye & $\kisc$ & $p$ & MLAE & canonical & channel & calls & Toffolis \\
    \midrule""")
for n in ALL:
    if n == CANDIDATES[0]:
        w(r"    \midrule")
    r = dc[n]
    d1, d2 = r.get("D1", {}), r.get("D2", {})
    ch = (d1.get("channel_route") or {}).get("toffoli")
    if d2.get("free"):
        c2, t2 = "--", "free"
    else:
        c2, t2 = sci(d2.get("calls_mlae")), sci(d2.get("toffoli_mlae"))
    w(f"    {lab(n)} & {sci(r['k_isc'])} & {sci(r['p'])} & {sci(d1.get('toffoli_mlae'))} & "
      f"{sci(d1.get('toffoli_canonical'))} & {sci(ch)} & {c2} & {t2} \\\\")
w(r"""    \bottomrule
  \end{tabular}
\end{table}
""")

# --- geometries -----------------------------------------------------------------------------------------------
w(r"\input{si_tn_shift}\clearpage  % written by paper/analysis/tn_shift_check.py")
w(r"\input{si_tda_full}\clearpage  % written by paper/analysis/tda_full_check.py")
w(r"\input{si_xc_conf}\clearpage  % written by paper/analysis/xc_conformer_table.py")
w(r"\input{si_fcwd}\clearpage  % written by paper/analysis/fcwd_table.py")
w(r"""\section{Cartesian coordinates}
GFN2-xTB ground-state geometries (\AA) used for every calculation in the paper (lowest conformer for the distyryl
dyes). The same files are in the repository (\texttt{validation/runs}, \texttt{screening/runs/candidates\_round1}).
\begin{multicols}{2}""")
for n in ALL:
    p = geo_path(n)
    lines = p.read_text().split("\n")
    nat = int(lines[0])
    w(f"\\paragraph{{{lab(n)}}} ({nat} atoms)\n{{\\scriptsize\\linespread{{1}}\\selectfont\\begin{{verbatim}}")
    for ln in lines[2:2 + nat]:
        s = ln.split()
        w(f"{s[0]:<3s}{float(s[1]):12.6f}{float(s[2]):12.6f}{float(s[3]):12.6f}")
    w("\\end{verbatim}}")
w("\\end{multicols}")

(ROOT / "paper" / "si_tables.tex").write_text("\n".join(out) + "\n")
print("written paper/si_tables.tex")
