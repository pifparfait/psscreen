"""Active-space convergence sweep: do S1-T1 and T2-S1 gaps approach experiment as the active space grows?

For one dye, runs a nested sequence of active spaces (default (8,8), (10,10), (12,12), (14,14), (16,16)). Each size
is warm-started from the orbitals of the previous one (the new active space adds the highest core and the lowest
virtual orbital to the previous one):
  * SA-CASSCF over the lowest singlets (default 3, equal weights) up to --casscf-max orbitals (default 14);
    larger sizes use CASCI on the orbitals of the previous size;
  * CASCI singlets and triplets (M_S = 1 sector, spin-fixed) on those orbitals;
  * SC-NEVPT2 for every root up to --nevpt2-max orbitals (default 12; the 4-RDM grows as n^8).
Results are written after every step to paper/references/results/NAME.cas_sweep.json, so a partial sweep is
usable. Natural-orbital occupations of the state-averaged density are stored for each size.

Usage (repository root): python paper/references/run_cas_sweep.py NAME [--sizes 8,10,12,14,16]
Environment: OMP_NUM_THREADS, PSS_MAX_MEMORY_MB (default 60000).
"""
import argparse
import json
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
import numpy as np  # noqa: E402
from pyscf import fci, gto, mcscf, mrpt, scf  # noqa: E402

from psscreen.backends import _ecp  # noqa: E402
from psscreen.geometry import read_xyz, to_pyscf_atom  # noqa: E402

EV = 27.211386245988
ap = argparse.ArgumentParser()
ap.add_argument("name")
ap.add_argument("--sizes", default="8,10,12,14,16", help="n for (n,n) active spaces, increasing")
ap.add_argument("--nsing", type=int, default=3)
ap.add_argument("--ntrip", type=int, default=3)
ap.add_argument("--casscf-max", type=int, default=14)
ap.add_argument("--nevpt2-max", type=int, default=12)
ap.add_argument("--max-macro", type=int, default=50)
a = ap.parse_args()
sizes = [int(x) for x in a.sizes.split(",")]
assert sizes == sorted(sizes)

geo = os.path.join(ROOT, "validation", "runs", f"{a.name}.xyz")
if not os.path.exists(geo):
    geo = os.path.join(ROOT, "validation", "runs", "set3", f"{a.name}_c0.xyz")
sym, xyz = read_xyz(geo)
atom = to_pyscf_atom(sym, xyz)
mem = float(os.environ.get("PSS_MAX_MEMORY_MB", 60000))
mol = gto.M(atom=atom, basis="def2-svp", ecp=_ecp(atom, "def2-svp"), verbose=4, max_memory=mem)
outdir = os.path.join(ROOT, "paper", "references", "results")
os.makedirs(outdir, exist_ok=True)
fout = os.path.join(outdir, f"{a.name}.cas_sweep.json")
out = {"name": a.name, "geometry": os.path.relpath(geo, ROOT), "basis": "def2-svp", "nsing": a.nsing,
       "ntrip": a.ntrip, "casscf_max": a.casscf_max, "nevpt2_max": a.nevpt2_max, "sizes": {}}


def save():
    json.dump(out, open(fout, "w"), indent=1)


def solver(spin, n):
    s = fci.direct_spin0.FCI(mol) if spin == 0 else fci.direct_spin1.FCI(mol)
    if n >= 16:
        s.max_space = 6           # keep the Davidson subspace small for ~1e8 determinants
    s.max_memory = mem
    return s


t0 = time.time()
mf = scf.RHF(mol).density_fit().run()
mo = mf.mo_coeff
out["E_RHF"] = float(mf.e_tot)
save()

for n in sizes:
    key = f"({n},{n})"
    rec = {"timings_s": {}}
    out["sizes"][key] = rec
    t = time.time()
    if n <= a.casscf_max:
        mc = mcscf.CASSCF(mf, n, n).density_fit()
        mc.fcisolver = solver(0, n)
        mc.fcisolver.nroots = a.nsing
        mc = mcscf.state_average_(mc, [1.0 / a.nsing] * a.nsing)
        mc.max_cycle_macro = a.max_macro
        mc.kernel(mo)
        rec["orbitals"] = f"SA-CASSCF({n},{n}) over {a.nsing} singlets"
        rec["sa_casscf_converged"] = bool(mc.converged)
        mo = mc.mo_coeff
    else:
        rec["orbitals"] = "orbitals of the previous size (CASCI only)"
    rec["timings_s"]["orbitals"] = time.time() - t
    save()

    t = time.time()
    cs = mcscf.CASCI(mf, n, n)
    cs.fcisolver = solver(0, n)
    cs.fcisolver.nroots = a.nsing
    cs.kernel(mo)
    ct = mcscf.CASCI(mf, n, n)
    ct.fcisolver = solver(2, n)
    ct.nelecas = ((n + 2) // 2, (n - 2) // 2)
    ct.fix_spin_(ss=2.0)
    ct.fcisolver.nroots = a.ntrip
    ct.kernel(mo)
    E0 = cs.e_tot[0]
    rec["E_CASCI_singlets_eV"] = [(e - E0) * EV for e in cs.e_tot]
    rec["E_CASCI_triplets_eV"] = [(e - E0) * EV for e in ct.e_tot]
    dm1 = sum(cs.fcisolver.make_rdm1(c, n, n) for c in cs.ci) / a.nsing
    rec["sa_natural_occupations"] = sorted(np.linalg.eigvalsh(dm1).tolist(), reverse=True)
    rec["timings_s"]["casci"] = time.time() - t
    print(a.name, key, "CASCI S:", [round(x, 3) for x in rec["E_CASCI_singlets_eV"]],
          "T:", [round(x, 3) for x in rec["E_CASCI_triplets_eV"]], flush=True)
    save()

    if n <= a.nevpt2_max:
        t = time.time()
        try:
            es = [cs.e_tot[i] + mrpt.NEVPT(cs, root=i).kernel() for i in range(a.nsing)]
            rec["E_NEVPT2_singlets_eV"] = [(e - es[0]) * EV for e in es]
            save()
            et = [ct.e_tot[i] + mrpt.NEVPT(ct, root=i).kernel() for i in range(a.ntrip)]
            rec["E_NEVPT2_triplets_eV"] = [(e - es[0]) * EV for e in et]
            print(a.name, key, "NEVPT2 S:", [round(x, 3) for x in rec["E_NEVPT2_singlets_eV"]],
                  "T:", [round(x, 3) for x in rec["E_NEVPT2_triplets_eV"]], flush=True)
        except Exception as exc:
            rec["nevpt2_error"] = repr(exc)
            print(a.name, key, "NEVPT2 failed:", exc, flush=True)
        rec["timings_s"]["nevpt2"] = time.time() - t
    save()
    print(a.name, key, "done", round(time.time() - t0), "s", flush=True)

out["timings_s_total"] = time.time() - t0
save()
