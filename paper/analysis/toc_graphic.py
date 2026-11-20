"""Table-of-contents graphic for JCTC (8.25 cm x 4.45 cm): the gap tolerance that a factor-2 ISC rate demands, and
the Toffoli cost of one circuit call versus one screening decision.
Reads paper/analysis/tolerance_curve.json and paper/quantum/results/decision_cost.json.
Usage (repository root): python paper/analysis/toc_graphic.py   -> paper/figures/toc.pdf, toc.png
"""
import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
tol = json.loads((ROOT / "paper" / "analysis" / "tolerance_curve.json").read_text())
dec = json.loads((ROOT / "paper" / "quantum" / "results" / "decision_cost.json").read_text())

plt.rcParams.update({"font.size": 6, "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
                     "xtick.major.size": 2, "ytick.major.size": 2, "font.family": "DejaVu Sans"})
fig, (a1, a2) = plt.subplots(1, 2, figsize=(8.25 / 2.54, 4.45 / 2.54), gridspec_kw={"width_ratios": [1.15, 1]})

# left: tolerance vs allowed error factor
f = np.array(tol["factors"])
curves = np.array([v["vs_factor"] for v in tol["dyes"].values()])
lo, hi, med = curves.min(0), curves.max(0), np.median(curves, 0)
a1.fill_between(f, lo, hi, color="#90cdf4", alpha=0.6, lw=0)
a1.plot(f, med, color="#2b6cb0", lw=1)
a1.set_xscale("log")
a1.set_xticks([1.25, 2, 5, 10])
a1.set_xticklabels(["1.25", "2", "5", "10"])
a1.minorticks_off()
a1.axvline(2, color="0.5", lw=0.5, ls=":")
a1.annotate("factor 2:\ngaps to ≈0.08 eV", xy=(2, 0.07), xytext=(2.6, 0.02), fontsize=5.5,
            arrowprops=dict(arrowstyle="-", lw=0.4, color="0.3"))
a1.set_xlabel(r"allowed error in $k_\mathrm{ISC}$ (factor)")
a1.set_ylabel("gap tolerance (eV)")
a1.set_ylim(0, 0.3)
a1.set_title("accuracy needed", fontsize=6.5, pad=2)

# right: Toffolis per call vs per decision
rows = [r for r in dec["rows"] if isinstance(r.get("D1"), dict) and "toffoli_mlae" in r["D1"]]
call = [r["D1"]["C_call"] for r in rows]
d2 = [r["D2"]["toffoli_mlae"] for r in rows if isinstance(r.get("D2"), dict) and "toffoli_mlae" in r["D2"]]
d1 = [r["D1"]["toffoli_mlae"] for r in rows]
rng = np.random.default_rng(1)
for x, ys, c in ((0, call, "#718096"), (1, d2, "#2f855a"), (2, d1, "#c05621")):
    a2.plot(x + 0.12 * rng.standard_normal(len(ys)), ys, "o", ms=2, color=c, alpha=0.8, mec="none")
a2.set_yscale("log")
a2.set_ylim(1e7, 1e20)
a2.set_yticks([1e8, 1e11, 1e14, 1e17, 1e20])
a2.minorticks_off()
a2.set_xticks([0, 1, 2])
a2.set_xticklabels(["one\ncall", r"$\Phi_T>0.5$?", r"$k_\mathrm{ISC}$" + "\nto 20%"])
a2.set_xlim(-0.5, 2.5)
a2.set_ylabel("Toffoli gates")
a2.set_title("cost of a decision", fontsize=6.5, pad=2)
for a in (a1, a2):
    a.spines[["top", "right"]].set_visible(False)
fig.tight_layout(pad=0.3, w_pad=0.8)
out = ROOT / "paper" / "figures"
fig.savefig(out / "toc.pdf")
fig.savefig(out / "toc.png", dpi=600)
print("written", out / "toc.pdf")
