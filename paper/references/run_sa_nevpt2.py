"""Multireference reference for S1, T1, T2, ...: SA-CASSCF orbitals + CASCI/SC-NEVPT2 singlets and triplets.

This is the classical bound on what an exact active-space solver (e.g. quantum phase estimation) plus a
perturbative dynamic-correlation correction would give. Compared with validation/run_casci_reference.py
it (i) optimizes the orbitals by state-averaged CASSCF over the lowest singlets instead of using RHF
orbitals, and (ii) also computes triplets, so the S1-T1 and S1-T2 gaps that control ISC can be checked.

Usage (from the repository root):
    python paper/references/run_sa_nevpt2.py NAME [--cas 12,11] [--nsing 3] [--ntrip 3] [--no-nevpt2]
NAME is a dye with a geometry in validation/runs/NAME.xyz or validation/runs/set3/NAME_c0.xyz.
Environment: OMP_NUM_THREADS, PSS_MAX_MEMORY_MB (default 20000).
Writes paper/references/results/NAME.sa_nevpt2_<ne>_<no>.json (updated after every step).
"""
import argparse
import json
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from pyscf import fci, gto, mcscf, mrpt, scf  # noqa: E402

from psscreen.backends import _ecp  # noqa: E402
from psscreen.geometry import read_xyz, to_pyscf_atom  # noqa: E402

EV = 27.211386245988

ap = argparse.ArgumentParser()
ap.add_argument("name")
ap.add_argument("--cas", default="12,11", help="electrons,orbitals")
ap.add_argument("--nsing", type=int, default=3, help="singlet roots (S0 included)")
ap.add_argument("--ntrip", type=int, default=3)
ap.add_argument("--no-nevpt2", action="store_true")
a = ap.parse_args()
ne, no = (int(x) for x in a.cas.split(","))

geo = os.path.join(ROOT, "validation", "runs", f"{a.name}.xyz")
if not os.path.exists(geo):
    geo = os.path.join(ROOT, "validation", "runs", "set3", f"{a.name}_c0.xyz")
sym, xyz = read_xyz(geo)
atom = to_pyscf_atom(sym, xyz)
mol = gto.M(atom=atom, basis="def2-svp", ecp=_ecp(atom, "def2-svp"), verbose=4,
            max_memory=float(os.environ.get("PSS_MAX_MEMORY_MB", 20000)))
outdir = os.path.join(ROOT, "paper", "references", "results")
os.makedirs(outdir, exist_ok=True)
fout = os.path.join(outdir, f"{a.name}.sa_nevpt2_{ne}_{no}.json")
out = {"name": a.name, "geometry": os.path.relpath(geo, ROOT), "cas": [ne, no], "basis": "def2-svp",
       "orbitals": f"SA-CASSCF({ne},{no}) over {a.nsing} singlets, density-fitted", "timings_s": {}}


def save():
    json.dump(out, open(fout, "w"), indent=1)


t0 = time.time()
mf = scf.RHF(mol).density_fit().run()
out["timings_s"]["rhf"] = time.time() - t0

# 1) state-averaged CASSCF orbitals (singlets only, equal weights)
t = time.time()
mc = mcscf.CASSCF(mf, no, ne).density_fit()
mc.fcisolver = fci.direct_spin0.FCI(mol)
mc.fcisolver.nroots = a.nsing
mc = mcscf.state_average_(mc, [1.0 / a.nsing] * a.nsing)
mc.max_cycle_macro = 100
mc.kernel()
out["sa_casscf_converged"] = bool(mc.converged)
out["timings_s"]["sa_casscf"] = time.time() - t
mo = mc.mo_coeff
save()


def casci(spin, nroots):
    c = mcscf.CASCI(mf, no, ne)
    if spin == 0:
        c.fcisolver = fci.direct_spin0.FCI(mol)
    else:
        c.fcisolver = fci.direct_spin1.FCI(mol)
        c.nelecas = ((ne + spin) // 2, (ne - spin) // 2)        # M_S = S, pure triplets in this sector
        c.fix_spin_(ss=spin / 2 * (spin / 2 + 1))
    c.fcisolver.nroots = nroots
    c.kernel(mo)
    return c


# 2) CASCI singlets and triplets on the SA-CASSCF orbitals
t = time.time()
cs = casci(0, a.nsing)
ct = casci(2, a.ntrip)
E0 = cs.e_tot[0]
out["E_CASCI_singlets_eV"] = [(e - E0) * EV for e in cs.e_tot]
out["E_CASCI_triplets_eV"] = [(e - E0) * EV for e in ct.e_tot]
out["timings_s"]["casci"] = time.time() - t
print(a.name, "CASCI S:", [round(x, 3) for x in out["E_CASCI_singlets_eV"]],
      "T:", [round(x, 3) for x in out["E_CASCI_triplets_eV"]], flush=True)
save()

# 3) SC-NEVPT2 for every root
if not a.no_nevpt2:
    t = time.time()
    try:
        es = [cs.e_tot[i] + mrpt.NEVPT(cs, root=i).kernel() for i in range(a.nsing)]
        out["E_NEVPT2_singlets_eV"] = [(e - es[0]) * EV for e in es]
        save()
        et = [ct.e_tot[i] + mrpt.NEVPT(ct, root=i).kernel() for i in range(a.ntrip)]
        out["E_NEVPT2_triplets_eV"] = [(e - es[0]) * EV for e in et]
        print(a.name, "NEVPT2 S:", [round(x, 3) for x in out["E_NEVPT2_singlets_eV"]],
              "T:", [round(x, 3) for x in out["E_NEVPT2_triplets_eV"]], flush=True)
    except Exception as exc:  # keep CASCI results if NEVPT2 fails (memory)
        out["nevpt2_error"] = repr(exc)
        print(a.name, "NEVPT2 failed:", exc, flush=True)
    out["timings_s"]["nevpt2"] = time.time() - t
out["timings_s"]["total"] = time.time() - t0
save()
print(a.name, "done", round(time.time() - t0), "s", flush=True)
