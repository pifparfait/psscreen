"""CLI de cribado. Ejemplo:
    python screen.py candidatos.csv --out resultados/
candidatos.csv con columnas: name,smiles
"""
import argparse, csv, json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from psscreen.pipeline import screen_molecule
from psscreen.backends import PySCFTDA

ap = argparse.ArgumentParser()
ap.add_argument("csv"); ap.add_argument("--out", default="runs")
ap.add_argument("--xc", default="b3lyp"); ap.add_argument("--basis", default="def2-svp")
ap.add_argument("--k-ic", type=float, default=1.5e8, help="conversión interna S1->S0 (s^-1), empírica")
ap.add_argument("--n-refr", type=float, default=1.496, help="índice de refracción del disolvente")
ap.add_argument("--window", default="650,850", help="ventana terapéutica en nm")
a = ap.parse_args()
win = tuple(float(v) for v in a.window.split(","))
rows = list(csv.DictReader(open(a.csv)))
res = []
for r in rows:
    o = screen_molecule(r["name"], r["smiles"], backend=PySCFTDA(xc=a.xc, basis=a.basis),
                        workdir=a.out, k_ic=a.k_ic, n_refr=a.n_refr, window=win)
    p = o["photophysics"]
    res.append((o["pdt_score"], r["name"], o["lambda_max_nm"], p["phi_ISC"], o["E_T1_eV"]))
    print(f"{r['name']}: λ={o['lambda_max_nm']:.0f} nm Φ_T={p['phi_ISC']:.2f} T1={o['E_T1_eV']:.2f} eV score={o['pdt_score']:.2f}", flush=True)
print("\nRanking PDT:")
for s, n, l, ph, t in sorted(res, reverse=True):
    print(f"  {n:20s} score={s:.2f}  λ={l:.0f} nm  Φ_T={ph:.2f}  T1={t:.2f} eV")
