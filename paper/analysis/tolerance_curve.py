"""Gap tolerance as a function of the allowed k_ISC factor and of the solvent reorganization energy.

Uses the same perturbation model as accuracy_budget.py (all triplets shifted rigidly relative to S1) and writes
paper/figures/tolerance_curve.pdf and paper/analysis/tolerance_curve.json.
Usage (from the repository root): python paper/analysis/tolerance_curve.py
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "paper" / "analysis"))
sys.path.insert(0, str(ROOT))
from accuracy_budget import CAND_DIR, CANDIDATES, DYES, load, tolerance  # noqa: E402
from psscreen.rates import MLJParams  # noqa: E402

FACTORS = [1.25, 1.5, 2, 3, 5, 10]
LAMBDAS = [0.05, 0.075, 0.10, 0.15, 0.20, 0.25]
SHOW = ["I2-BDP", "I3-BDP", "I2-BDP-0", "BDP", "azaBDP_Ph4", "azaBDP_Br2"]


def main():
    data = {}
    for n in DYES + CANDIDATES:
        d = load(n) or load(n, CAND_DIR)
        data[n] = {"vs_factor": [tolerance(d, f, "all", MLJParams(), dmax=0.8, n=801) for f in FACTORS],
                   "vs_lambda": [tolerance(d, 2.0, "all", MLJParams(lambda_s=l), dmax=0.8, n=801)
                                 for l in LAMBDAS]}
        print(n, data[n], flush=True)
    (ROOT / "paper" / "analysis" / "tolerance_curve.json").write_text(
        json.dumps({"factors": FACTORS, "lambda_s": LAMBDAS, "dyes": data}, indent=1))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 2, figsize=(6.8, 2.8))
    val = [n for n in DYES]
    for ax, key, xs, xl in ((axs[0], "vs_factor", FACTORS, r"allowed error factor in $k_\mathrm{ISC}$"),
                            (axs[1], "vs_lambda", LAMBDAS, r"$\lambda_s$ (eV)")):
        arr = np.array([data[n][key] for n in val])
        ax.fill_between(xs, arr.min(0), arr.max(0), color="0.85", label="nine validation dyes (range)")
        for n, c in zip(SHOW, ["#2b6cb0", "#2c5282", "#c05621", "0.4", "#2f855a", "#9b2c2c"]):
            ax.plot(xs, data[n][key], "o-", ms=2.5, lw=1, color=c, label=n.replace("_", " "))
        ax.set_xlabel(xl)
        ax.set_ylabel("tolerated S$_1$–T gap error (eV)")
        if key == "vs_factor":
            ax.set_xscale("log")
            from matplotlib.ticker import NullFormatter, NullLocator
            ax.xaxis.set_minor_locator(NullLocator())
            ax.xaxis.set_minor_formatter(NullFormatter())
            ax.set_xticks([1.5, 2, 3, 5, 10])
            ax.set_xticklabels(["1.5", "2", "3", "5", "10"])
            ax.axvline(2, color="0.5", ls=":", lw=0.8)
        else:
            ax.axvline(0.10, color="0.5", ls=":", lw=0.8)
    axs[0].legend(frameon=False, fontsize=6)
    axs[0].set_title("(a) $\\lambda_s$ = 0.10 eV", fontsize=8)
    axs[1].set_title("(b) factor 2", fontsize=8)
    fig.tight_layout()
    fig.savefig(ROOT / "paper" / "figures" / "tolerance_curve.pdf")


if __name__ == "__main__":
    main()
