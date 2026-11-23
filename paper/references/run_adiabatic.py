"""Vertical-to-adiabatic corrections for the S1 and T1 energies, with the same DeltaSCF treatment as the pipeline.

The pipeline and all multireference results of the paper use vertical energies at the GFN2-xTB ground-state geometry,
whereas the MLJ rate expression assumes adiabatic gaps. This script optimizes, starting from the xTB geometry,
  S0 : RKS B3LYP,
  T1 : UKS B3LYP (M_S = 1, HOMO->LUMO occupation held by MOM),
  S1 : DeltaSCF sum rule E_S1 = 2 E_mix - E_T1, with E_mix the M_S = 0 HOMO->LUMO mixed determinant (MOM); its
       gradient is 2 g_mix - g_T1,
with ASE's BFGS optimizer driving PySCF energies and analytic gradients (B3LYP/def2-SVP, density fitting, def2 ECP on
iodine, grid level 2, as in psscreen.deltascf). This is the treatment that gives the pipeline's S1 and T1, and it avoids
TDDFT gradients, which PySCF does not provide with density fitting. At every minimum it evaluates the S0, T1 and S1
energies and reports the relaxation energies
    lambda_S1 = E_S1(S0 min) - E_S1(S1 min),   lambda_T1 = E_T1(S0 min) - E_T1(T1 min),
so that the adiabatic S1-T1 gap is the vertical gap at the S0 minimum minus lambda_S1 plus lambda_T1, plus the effect of
replacing the xTB geometry by the DFT S0 minimum. Each optimization is check-pointed: the last geometry is written to
paper/references/results/NAME.adiabatic_<state>.xyz after every step and reused on restart.

Usage (repository root): PYTHONUNBUFFERED=1 python paper/references/run_adiabatic.py NAME [--states S0 T1 S1]
Environment: OMP_NUM_THREADS, PSS_MAX_MEMORY_MB (default 40000).
Writes paper/references/results/NAME.adiabatic.json.
"""
import argparse
import json
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402
from ase import Atoms  # noqa: E402
from ase.calculators.calculator import Calculator, all_changes  # noqa: E402
from ase.optimize import BFGS  # noqa: E402
from pyscf import dft, gto, scf  # noqa: E402

from psscreen.backends import _ecp  # noqa: E402
from psscreen.geometry import read_xyz, to_pyscf_atom, write_xyz  # noqa: E402

HARTREE_EV = 27.211386245988
BOHR_A = 0.52917721092

ap = argparse.ArgumentParser()
ap.add_argument("name")
ap.add_argument("--states", nargs="+", default=["S0", "T1", "S1"])
ap.add_argument("--basis", default="def2-svp")
ap.add_argument("--xc", default="b3lyp")
ap.add_argument("--fmax", type=float, default=0.03, help="eV/Angstrom")
ap.add_argument("--max-steps", type=int, default=200)
a = ap.parse_args()

RES = os.path.join(ROOT, "paper", "references", "results")
OUT = os.path.join(RES, f"{a.name}.adiabatic.json")
MEM = int(os.environ.get("PSS_MAX_MEMORY_MB", 40000))
res = json.load(open(OUT)) if os.path.exists(OUT) else {"name": a.name, "method": f"{a.xc}/{a.basis}, DF; T1 UKS (MOM), S1 DeltaSCF sum rule",
                                                         "minima": {}, "timings_s": {}}


def save():
    os.makedirs(RES, exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(res, f, indent=1)


def build(sym, xyz, spin=0):
    atom = to_pyscf_atom(sym, xyz)
    mol = gto.M(atom=atom, basis=a.basis, ecp=_ecp(atom, a.basis), spin=spin, verbose=0, max_memory=MEM)
    mf = (dft.UKS if spin else dft.RKS)(mol, xc=a.xc).density_fit()
    mf.conv_tol = 1e-9
    mf.max_cycle = 150
    mf.grids.level = 2
    return mol, mf


_prev = {}      # converged orbitals per state, reused as guess and MOM reference at the next geometry


def _occupations(nmo, nocc, kind):
    oa, ob = np.zeros(nmo), np.zeros(nmo)
    if kind == "T":
        oa[:nocc + 1] = 1
        ob[:nocc - 1] = 1
    else:                          # mixed determinant: alpha HOMO -> alpha LUMO
        oa[:nocc - 1] = 1
        oa[nocc] = 1
        ob[:nocc] = 1
    return oa, ob


def _uks_state(sym, xyz, kind):
    """UKS with MOM for the triplet (kind 'T') or the mixed determinant (kind 'M'); returns the converged object."""
    mol, mf = build(sym, xyz, spin=2 if kind == "T" else 0)
    if kind == "M":
        mf = dft.UKS(mol, xc=a.xc).density_fit()
        mf.conv_tol, mf.max_cycle, mf.grids.level = 1e-9, 150, 2
    if kind in _prev:
        mo0, occ = _prev[kind]
    else:
        _, m0 = build(sym, xyz, 0)
        m0.kernel()
        nocc = mol.nelectron // 2
        mo0 = (m0.mo_coeff, m0.mo_coeff)
        occ = _occupations(len(m0.mo_energy), nocc, kind)
    mf = scf.addons.mom_occ(mf, mo0, occ)
    dm0 = np.array([(mo0[s] * occ[s]) @ mo0[s].T for s in (0, 1)])
    mf.kernel(dm0)
    if not mf.converged:
        # second attempt: level shift and damping from the last density, then a polish without level shift
        mol2, mf2 = build(sym, xyz, spin=2 if kind == "T" else 0)
        if kind == "M":
            mf2 = dft.UKS(mol2, xc=a.xc).density_fit()
            mf2.grids.level = 2
        mf2.conv_tol, mf2.max_cycle = 1e-9, 400
        mf2.level_shift, mf2.damp, mf2.diis_start_cycle = 0.25, 0.5, 8
        mf2 = scf.addons.mom_occ(mf2, mo0, occ)
        mf2.kernel(mf.make_rdm1())
        if mf2.converged:
            mol3, mf3 = build(sym, xyz, spin=2 if kind == "T" else 0)
            if kind == "M":
                mf3 = dft.UKS(mol3, xc=a.xc).density_fit()
                mf3.grids.level = 2
            mf3.conv_tol, mf3.max_cycle = 1e-9, 100
            mf3 = scf.addons.mom_occ(mf3, mf2.mo_coeff, mf2.mo_occ)
            mf3.kernel(mf2.make_rdm1())
            mf = mf3 if mf3.converged else mf2
    assert mf.converged, f"UKS ({kind}) not converged"
    _prev[kind] = ((mf.mo_coeff[0], mf.mo_coeff[1]), (mf.mo_occ[0], mf.mo_occ[1]))
    return mf


def energy_grad(sym, xyz, state, grad=True):
    """Total energy (Ha) and gradient (Ha/Bohr) of S0, T1 or S1 (DeltaSCF sum rule)."""
    if state == "S0":
        mol, mf = build(sym, xyz, 0)
        dm0 = _prev.get("S0")
        mf.kernel(dm0=dm0)
        assert mf.converged, "RKS not converged"
        _prev["S0"] = mf.make_rdm1()
        return mf.e_tot, (mf.nuc_grad_method().kernel() if grad else None)
    mt = _uks_state(sym, xyz, "T")
    if state == "T1":
        return mt.e_tot, (mt.nuc_grad_method().kernel() if grad else None)
    mm = _uks_state(sym, xyz, "M")
    e = 2 * mm.e_tot - mt.e_tot
    if not grad:
        return e, None
    return e, 2 * mm.nuc_grad_method().kernel() - mt.nuc_grad_method().kernel()


class PySCFCalc(Calculator):
    implemented_properties = ["energy", "forces"]

    def __init__(self, state, sym, ckpt):
        super().__init__()
        self.state, self.sym, self.ckpt = state, sym, ckpt

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        xyz = self.atoms.get_positions()
        e, g = energy_grad(self.sym, xyz, self.state)
        self.results = {"energy": e * HARTREE_EV, "forces": -np.asarray(g) * HARTREE_EV / BOHR_A}
        write_xyz(self.ckpt, self.sym, xyz, f"{a.name} {self.state} E={e:.10f}")


def all_energies(sym, xyz):
    out = {}
    for st in ("S0", "T1", "S1"):
        out[f"E_{st}_Ha"] = energy_grad(sym, xyz, st, grad=False)[0]
    out["S1_minus_T1_eV"] = (out["E_S1_Ha"] - out["E_T1_Ha"]) * HARTREE_EV
    out["S1_vert_eV"] = (out["E_S1_Ha"] - out["E_S0_Ha"]) * HARTREE_EV
    out["T1_vert_eV"] = (out["E_T1_Ha"] - out["E_S0_Ha"]) * HARTREE_EV
    return out


t0 = time.time()
geo = os.path.join(ROOT, "validation", "runs", f"{a.name}.xyz")
for alt in (os.path.join(ROOT, "validation", "runs", "set3", f"{a.name}_c0.xyz"),
            os.path.join(ROOT, "screening", "runs", "candidates_round1", f"{a.name}.xyz")):
    if not os.path.exists(geo):
        geo = alt
sym, xyz0 = read_xyz(geo)
sym = [s.capitalize() for s in sym]
xyz0 = np.asarray(xyz0, float)

if "xtb" not in res["minima"]:
    res["minima"]["xtb"] = all_energies(sym, xyz0)
    save()
print("xTB geometry:", res["minima"]["xtb"], flush=True)

for st in a.states:
    if st in res["minima"]:
        print(st, "already done", flush=True)
        continue
    ckpt = os.path.join(RES, f"{a.name}.adiabatic_{st}.xyz")
    if os.path.exists(ckpt):
        start = np.asarray(read_xyz(ckpt)[1], float)
    elif st != "S0" and os.path.exists(os.path.join(RES, f"{a.name}.adiabatic_S0.xyz")) and "S0" in res["minima"]:
        start = np.asarray(read_xyz(os.path.join(RES, f"{a.name}.adiabatic_S0.xyz"))[1], float)
    else:
        start = xyz0
    atoms = Atoms(symbols=sym, positions=start)
    atoms.calc = PySCFCalc(st, sym, ckpt)
    t1 = time.time()
    opt = BFGS(atoms, logfile=os.path.join(RES, f"{a.name}.adiabatic_{st}.log"))
    conv = opt.run(fmax=a.fmax, steps=a.max_steps)
    xyz = atoms.get_positions()
    write_xyz(ckpt, sym, xyz, f"{a.name} {st} minimum")
    _prev.clear()
    res["minima"][st] = all_energies(sym, xyz)
    res["minima"][st].update({"converged": bool(conv), "steps": opt.nsteps,
                              "rmsd_from_xtb_A": float(np.sqrt(((xyz - xyz0) ** 2).sum(1).mean()))})
    res["timings_s"][st] = time.time() - t1
    save()
    print(st, res["minima"][st], flush=True)

m = res["minima"]
if all(k in m for k in ("S0", "T1", "S1")):
    lam_s = (m["S0"]["E_S1_Ha"] - m["S1"]["E_S1_Ha"]) * HARTREE_EV
    lam_t = (m["S0"]["E_T1_Ha"] - m["T1"]["E_T1_Ha"]) * HARTREE_EV
    res["summary"] = {
        "lambda_S1_eV": lam_s, "lambda_T1_eV": lam_t,
        "gap_vertical_xtb_eV": m["xtb"]["S1_minus_T1_eV"],
        "gap_vertical_dft_S0_eV": m["S0"]["S1_minus_T1_eV"],
        "gap_adiabatic_eV": (m["S1"]["E_S1_Ha"] - m["T1"]["E_T1_Ha"]) * HARTREE_EV,
        "correction_vertical_to_adiabatic_eV": lam_t - lam_s,
        "S1_vertical_minus_adiabatic_eV": lam_s,
    }
    print("summary:", res["summary"], flush=True)
res["timings_s"]["total"] = res["timings_s"].get("total", 0) + time.time() - t0
save()
print(a.name, "done", round(time.time() - t0), "s")
