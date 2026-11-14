"""Generate the figures of the paper draft from the validation JSON files.

Usage (from the repository root):  python paper/make_figures.py
Writes PDF figures to paper/figures/. Molecules without a result file are skipped,
so the styryl dyes appear automatically once their (trans) results are in validation/runs/.
"""
from __future__ import annotations
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "validation" / "runs"
OUT = ROOT / "paper" / "figures"
EV_NM = 1239.84198

REF = json.loads((ROOT / "validation" / "references.json").read_text())
ABS_SET = ["BODIPY_core", "BDP", "I-BDP", "I2-BDP", "I3-BDP", "azaBDP_Ph4",
           "BDP-0", "I-BDP-0", "I2-BDP-0"]
LABEL = {"BODIPY_core": "BODIPY core", "azaBDP_Ph4": "aza-BODIPY-Ph$_4$"}


def load(name):
    p = RUNS / f"{name}.json"
    return json.loads(p.read_text()) if p.exists() else None


def fig_absorption():
    xs, ys, names = [], [], []
    for n in ABS_SET:
        d = load(n)
        if d is None:
            continue
        xs.append(EV_NM / REF[n]["lambda_abs_nm"])
        ys.append(d["E_S1_eV"])
        names.append(LABEL.get(n, n))
    fig, ax = plt.subplots(figsize=(3.4, 3.2))
    lo, hi = 1.7, 2.6
    ax.plot([lo, hi], [lo, hi], color="0.6", lw=0.8, ls="--", label="ideal")
    ax.scatter(xs, ys, s=22, color="#2b6cb0", zorder=3, label=r"TDDFT+$\Delta$SCF")
    offs = {"BODIPY core": (-10, 6), "BDP": (6, 2), "I-BDP": (6, -8), "I2-BDP": (-36, 4),
            "I3-BDP": (6, -8)}
    for x, y, t in zip(xs, ys, names):
        ax.annotate(t, (x, y), textcoords="offset points", xytext=offs.get(t, (5, -9)), fontsize=6.5)
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_xlabel(r"$E_\mathrm{abs}$ experiment (eV)")
    ax.set_ylabel(r"$E(\mathrm{S}_1)$ calculated (eV)")
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    fig.tight_layout()
    fig.savefig(OUT / "absorption_parity.pdf")


def fig_cas_vs_exp():
    names = ["BODIPY_core", "BDP", "I2-BDP", "azaBDP_Ph4"]
    rows = []
    for n in names:
        c = json.loads((RUNS / f"{n}.casci_ref.json").read_text())
        d = load(n)
        rows.append((LABEL.get(n, n), EV_NM / REF[n]["lambda_abs_nm"], c["E_CASCI_eV"][1],
                     c["E_NEVPT2_eV"][1], d["E_S1_eV"]))
    fig, ax = plt.subplots(figsize=(3.4, 3.0))
    x = range(len(rows))
    w = 0.25
    for k, (lab, col) in enumerate([("CASCI(12,11)", "#c05621"), ("SC-NEVPT2", "#6b46c1"),
                                     (r"TDDFT+$\Delta$SCF", "#2b6cb0")]):
        ax.bar([i + (k - 1) * w for i in x], [r[2 + k] - r[1] for r in rows], w, label=lab, color=col)
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xticks(list(x)); ax.set_xticklabels([r[0] for r in rows], fontsize=7, rotation=15)
    ax.set_ylabel(r"$E(\mathrm{S}_1)$ calc $-$ exp (eV)")
    ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "active_space_error.pdf")


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    fig_absorption()
    fig_cas_vs_exp()
    print("figures written to", OUT)
