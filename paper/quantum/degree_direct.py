"""Direct filter degrees down to the polynomial precisions that decision_cost.py actually uses.

degree_check.py fitted the degree of the Chebyshev filter sqrt(rho_MLJ/rho_max) at eps = 1e-2, 2e-3 and 1e-4 only,
and decision_cost.py extrapolated a linear law in ln(1/eps) down to eps ~ 1e-9. Here the smallest even degree that
reaches a uniform error eps is found directly for eps = 1e-2 ... 1e-7 at alpha/Delta = 200 and 400 (the degree is
linear in alpha/Delta, degree_check.py), with Chebyshev coefficients from a DCT so that degrees of 1e4 are cheap.
Below eps ~ 1e-8 the target sqrt(rho) is limited by double-precision noise in rho (1e-16 relative, so 1e-8 after the
square root) and the degree is not meaningful. decision_cost.py interpolates the slope d/(alpha/Delta) in ln(1/eps)
over this grid and extrapolates linearly from the last two points below 1e-7 (needed only for D1 of dyes with
k_ISC < 1e6 s^-1, which are not reported).

Usage (from the repository root): python paper/quantum/degree_direct.py
Writes paper/quantum/results/degree_direct.json.
"""
import json
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "paper", "quantum"))
import numpy as np  # noqa: E402
from numpy.polynomial import chebyshev as C  # noqa: E402
from scipy.fft import dct  # noqa: E402

from degree_check import DELTA_HA, GAP_HA, f_factory  # noqa: E402


def cheb_coeffs(f, deg):
    n = deg + 1
    x = np.cos(np.pi * (np.arange(n) + 0.5) / n)
    c = dct(f(x), type=2) / n
    c[0] /= 2
    return c


def err_at(f, deg, xs, fx):
    return float(np.max(np.abs(C.chebval(xs, cheb_coeffs(f, deg)) - fx)))


def min_degree(alpha, eps):
    f = f_factory(alpha)
    x0 = GAP_HA / alpha
    w = 40 * DELTA_HA / alpha
    xs = np.clip(np.concatenate([np.linspace(-1, 1, 20001), np.linspace(x0 - w, x0 + 3 * w, 40001)]), -1, 1)
    fx = f(xs)
    lo, hi = 2, 64
    while err_at(f, hi, xs, fx) > eps:
        lo, hi = hi, 2 * hi
    while hi - lo > 2:
        mid = (lo + hi) // 2
        mid += mid % 2
        if err_at(f, mid, xs, fx) <= eps:
            hi = mid
        else:
            lo = mid
    return hi


def main():
    epss = [1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6, 1e-6, 3e-7, 1e-7]
    out = {"eps": epss, "ratios": {}}
    for r in (200, 400):
        ds = []
        for eps in epss:
            d = min_degree(r * DELTA_HA, eps)
            ds.append(d)
            print(f"alpha/Delta={r} eps={eps:.0e} degree={d} slope={d/r:.3f}", flush=True)
        out["ratios"][str(r)] = {"degree": ds, "slope": [d / r for d in ds]}
    json.dump(out, open(os.path.join(ROOT, "paper", "quantum", "results", "degree_direct.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
