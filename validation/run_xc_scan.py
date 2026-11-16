"""Functional scan on the validation dyes, reusing the GFN2-xTB geometries already in runs/.

Usage: python run_xc_scan.py XC [names...]     e.g.  python run_xc_scan.py bhandhlyp
Results go to runs_xc/<xc>/<name>.json. The TDDFT code supports global hybrids only
(B3LYP 20 % HF, PBE0 25 %, BHandHLYP 50 %), so the HF-exchange fraction is scanned
instead of using a range-separated functional.
"""
import sys, os, json
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(here)); sys.path.insert(0, here)
from psscreen.pipeline import screen_molecule
from psscreen.backends import PySCFTDA

DEFAULT = ["BODIPY_core", "BDP", "I-BDP", "I2-BDP", "I3-BDP", "azaBDP_Ph4",
           "BDP-0", "I-BDP-0", "I2-BDP-0"]
xc = sys.argv[1]
names = sys.argv[2:] or DEFAULT
out = os.path.join(here, "runs_xc", xc)
os.makedirs(out, exist_ok=True)
for n in names:
    if os.path.exists(os.path.join(out, f"{n}.json")):
        continue
    xyz = os.path.join(here, "runs", f"{n}.xyz")
    if not os.path.exists(xyz):  # flexible dyes: lowest-energy conformer
        xyz = os.path.join(here, "runs", "set3", f"{n}_c0.xyz")
    be = PySCFTDA(xc=xc, basis="def2-svp", nstates=5,
                  ewin_singlets=float(os.environ.get("PSS_EWIN_SINGLETS", 60)),
                  max_memory=int(os.environ.get("PSS_MAX_MEMORY_MB", 100000)))
    try:
        r = screen_molecule(n, xyz_path=xyz, backend=be, workdir=out, n_refr=1.496)
    except Exception as exc:
        print(f"{n} [{xc}] FAILED: {exc}", flush=True); continue
    p = r["photophysics"]
    print(f"{n} [{xc}]: lam={r['lambda_max_nm']:.0f} nm S1={r['E_S1_eV']:.2f} "
          f"kISC={p['k_isc']:.2e} phiT={p['phi_ISC']:.2f} t={r['wall_s']['total']:.0f}s", flush=True)
