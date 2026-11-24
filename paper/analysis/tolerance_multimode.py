"""Factor-2 gap tolerance with the multimode line shape (round-4 review).

For the S1<->T1 Huang-Rhys factors of BDP and I2-BDP (S0 normal modes, fcwd_multimode.py) the script
computes the Kubo density of states rho(dE) and, at each gap dE, the smallest rigid shift |delta| for which
rho(dE + delta)/rho(dE) leaves [1/2, 2]; the same for the single-mode MLJ model of the pipeline. Evaluated
at the pipeline T1 gaps and at the near-resonant uphill gap of aza-BDP-Br2 (+0.107 eV, using the I2-BDP
modes as a proxy, since no aza-BDP-Br2 Hessian is available).
Usage: python paper/analysis/tolerance_multimode.py  ->  paper/analysis/tolerance_multimode.{json,md}
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import numpy as np  # noqa: E402
from fcwd_multimode import (CM1_HA, RES, ROOT, huang_rhys, kabsch, normal_modes, read_xyz,  # noqa: E402
                            rho_multimode, MASS, BOHR_A)
from psscreen.rates import MLJParams, fcwd_mlj  # noqa: E402

GRID = np.round(np.arange(-1.2, 0.6001, 0.002), 6)


def s1t1_huang_rhys(name):
    hz = np.load(os.path.join(RES, f"{name}.s0_hessian.npz"))
    sym = [str(s) for s in hz["symbols"]]
    x0 = np.asarray(hz["xyz_A"], float)
    omega, L, sm, P = normal_modes(hz["H"], x0, sym)
    m_at = np.array([MASS[s] for s in sym])
    q = {}
    for st in ("S1", "T1"):
        xs = kabsch(x0, np.asarray(read_xyz(os.path.join(RES, f"{name}.adiabatic_{st}.xyz"))[1], float), m_at)
        q[st] = L.T @ (P @ ((xs - x0).ravel() / BOHR_A * sm))
    S, _ = huang_rhys(omega, q["S1"] - q["T1"])
    return omega, S


def tolerance(rho, dE):
    r0 = np.interp(dE, GRID, rho)
    best = np.inf
    for sgn in (1, -1):
        for d in np.arange(0.001, 0.6, 0.001):
            r = np.interp(dE + sgn * d, GRID, rho)
            if not (0.5 <= r / r0 <= 2.0):
                best = min(best, d)
                break
    return float(best)


def main():
    rows = []
    mlj = {ls: np.array([fcwd_mlj(float(x), MLJParams(lambda_s=ls)) for x in GRID]) for ls in (0.05, 0.10, 0.20)}
    targets = [("BDP", "BDP", -0.718), ("I2-BDP", "I2-BDP", -0.608), ("aza-BDP-Br2 (I2-BDP modes)", "I2-BDP", 0.107)]
    cache = {}
    for label, modes, dE in targets:
        if modes not in cache:
            omega, S = s1t1_huang_rhys(modes)
            cache[modes] = {ls: rho_multimode(GRID, omega, S, ls)[0] for ls in (0.0, 0.05, 0.10)}
        row = {"case": label, "dE_eV": dE,
               "multimode": {f"{ls:.2f}": tolerance(r, dE) for ls, r in cache[modes].items()},
               "mlj": {f"{ls:.2f}": tolerance(r, dE) for ls, r in mlj.items()}}
        rows.append(row)
    out = os.path.join(ROOT, "paper", "analysis", "tolerance_multimode")
    json.dump(rows, open(out + ".json", "w"), indent=1)
    L = ["# Factor-2 gap tolerance (eV): multimode Kubo vs single-mode MLJ", "",
         "| case | dE | multimode lam_out 0 / 0.05 / 0.10 | MLJ lam_s 0.05 / 0.10 / 0.20 |", "|---|---|---|---|"]
    for r in rows:
        m, s = r["multimode"], r["mlj"]
        L.append(f"| {r['case']} | {r['dE_eV']:+.3f} | {m['0.00']:.3f} / {m['0.05']:.3f} / {m['0.10']:.3f} | "
                 f"{s['0.05']:.3f} / {s['0.10']:.3f} / {s['0.20']:.3f} |")
    open(out + ".md", "w").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
