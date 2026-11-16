"""Compare functionals: absorption error, kISC and Phi_T vs experiment. Reads runs/ (B3LYP) and runs_xc/<xc>/."""
import json, os, glob
here = os.path.dirname(os.path.abspath(__file__))
ref = json.load(open(os.path.join(here, "references.json")))
EV = 1239.84198
dirs = {"b3lyp": None}
for d in sorted(glob.glob(os.path.join(here, "runs_xc", "*"))):
    if not os.path.isdir(d): continue
    dirs[os.path.basename(d)] = d
def load(xc, n):
    if xc == "b3lyp":
        for c in (f"{n}.json", os.path.join("set3", f"{n}_c0.json")):
            p = os.path.join(here, "runs", c)
            if os.path.exists(p): return json.load(open(p))
        return None
    p = os.path.join(dirs[xc], f"{n}.json")
    return json.load(open(p)) if os.path.exists(p) else None
print("| Dye | exp λ (nm) / kISC (s⁻¹) | " + " | ".join(f"{x}: ΔE (eV) / kISC / Φ_T" for x in dirs) + " |")
print("|---|---|" + "---|" * len(dirs))
for n, e in ref.items():
    cells = []
    for xc in dirs:
        r = load(xc, n)
        if r is None: cells.append("—"); continue
        p = r["photophysics"]
        cells.append(f"{r['E_S1_eV'] - EV / e['lambda_abs_nm']:+.2f} / {p['k_isc']:.1e} / {p['phi_ISC']:.2f}")
    k = e.get("k_isc_s"); k = f"{k:.1e}" if isinstance(k, float) else "—"
    print(f"| {n} | {e['lambda_abs_nm']} / {k} | " + " | ".join(cells) + " |")
