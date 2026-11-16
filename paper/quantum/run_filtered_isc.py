"""Gap-resolved ISC observable (psscreen/isc_filtered.py) on validation dyes, simulated with PennyLane.

Usage (from the repository root):
    python paper/quantum/run_filtered_isc.py NAME [NCAS ...] [--no-quantum]
NAME is a validation dye with validation/runs/NAME.xyz, or H2CO for a quick test.
Writes paper/quantum/results/NAME.filtered_isc.json.
"""
import json
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
import numpy as np  # noqa: E402
from pyscf import dft, gto  # noqa: E402

from psscreen.backends import _ecp  # noqa: E402
from psscreen.geometry import read_xyz, to_pyscf_atom  # noqa: E402
from psscreen.isc_filtered import filtered_isc  # noqa: E402
from psscreen.quantum_pl import active_space, eigenstates, qubit_hamiltonian  # noqa: E402
from psscreen.rates import MLJParams  # noqa: E402

H2CO = "C 0 0 0; O 0 0 1.208; H 0 0.943 -0.587; H 0 -0.943 -0.587"

args = [a for a in sys.argv[1:] if not a.startswith("--")]
run_q = "--no-quantum" not in sys.argv
name = args[0]
ncas_list = [int(v) for v in args[1:]] or [4]
if name == "H2CO":
    atom = H2CO
else:
    sym, xyz = read_xyz(os.path.join(ROOT, "validation", "runs", f"{name}.xyz"))
    atom = to_pyscf_atom(sym, xyz)
basis = "def2-svp"
mol = gto.M(atom=atom, basis=basis, ecp=_ecp(atom, basis), verbose=0)
t0 = time.time()
mf = dft.RKS(mol, xc="b3lyp").density_fit()
mf.grids.level = 2
mf.kernel()
print(name, f"SCF {time.time()-t0:.0f}s", flush=True)
outdir = os.path.join(ROOT, "paper", "quantum", "results")
os.makedirs(outdir, exist_ok=True)
out = {"name": name, "cas": {}}
for n in ncas_list:
    t1 = time.time()
    asp = active_space(mf, n, n, basis)
    st = eigenstates(asp, qubit_hamiltonian(asp))
    res = {}
    for lab, p in {"lambda_s=0.10": MLJParams(), "lambda_s=0.20": MLJParams(lambda_s=0.20)}.items():
        r = filtered_isc(asp, st, p, run_quantum=run_q and lab == "lambda_s=0.10")
        res[lab] = r
        print(name, n, lab, f"deg={r['degree']} err={r['poly_err']:.1e} k_qsvt={r['k_isc_qsvt']} "
              f"k_poly={r['k_isc_poly']:.3e} k_exact={r['k_isc_exact']:.3e} proxy_T1_V={r['proxy_T1_V_cm']:.3f} cm-1 "
              f"dE_T1={r['dE_T1_eV']:+.2f} eV", flush=True)
    res["E_S_eV"] = st["E_S"][:4].tolist()
    res["E_T_eV"] = st["E_T"][:4].tolist()
    res["wall_s"] = time.time() - t1
    out["cas"][n] = res
    json.dump(out, open(os.path.join(outdir, f"{name}.filtered_isc.json"), "w"), indent=1, default=float)
