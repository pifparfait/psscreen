"""All-electron spin-orbit benchmark of the effective SOC operator used by psscreen.

psscreen computes SOC with an effective one-electron operator whose per-element scale reproduces the atomic valence
spin-orbit constant (psscreen/soc.py); for iodine it acts on ECP pseudo-orbitals. This script recomputes the same
S1-T_n couplings with a standard treatment:
  * all-electron scalar-relativistic reference: spin-free X2C-1e B3LYP, def2-SVP on light atoms and an uncontracted
    DZP-DKH basis on Br/I (default; --heavy-basis to change it);
  * TDA singlets and triplets (PySCF tdscf);
  * Breit-Pauli spin-orbit mean-field operator (SOMF): one-electron term plus the Coulomb and exchange mean-field
    two-electron terms built from the ground-state density,
        h_SO = (alpha^2/4) [ h1(pnucxp) + J(p1vxp1) - 3/2 (K + K') ] ,
    without picture-change correction. The one-electron-only (BP1e) couplings are reported as well.
The matrix elements reuse psscreen.soc (spin-adapted TDA states, sum over M). The script first checks the operator on
the free atoms (valence np spin-orbit constant versus experiment), then runs the molecule.

Usage (repository root): python paper/references/run_soc_benchmark.py NAME [--nsing 2] [--ntrip 5]
Environment: OMP_NUM_THREADS, PSS_MAX_MEMORY_MB (default 40000).
Writes paper/references/results/NAME.soc_benchmark.json.
"""
import argparse
import json
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402
from pyscf import dft, gto, scf, tdscf  # noqa: E402
from pyscf.scf import jk  # noqa: E402

from psscreen.geometry import read_xyz, to_pyscf_atom  # noqa: E402
from psscreen.soc import CM1_PER_HARTREE, ZETA_EXP_CM, _me, _spin_blocks, _state  # noqa: E402

ALPHA = 1 / 137.035999084
EV = 27.211386245988

ap = argparse.ArgumentParser()
ap.add_argument("name")
ap.add_argument("--nsing", type=int, default=2)
ap.add_argument("--ntrip", type=int, default=5)
ap.add_argument("--light-basis", default="def2-svp")
ap.add_argument("--heavy-basis", default="unc-dzp-dkh", help="basis for Br/I; prefix unc- to decontract")
ap.add_argument("--xc", default="b3lyp")
ap.add_argument("--atom-check-only", action="store_true")
ap.add_argument("--effective", action="store_true",
                help="also evaluate the effective one-electron operator of psscreen (per-element 1/r_A terms scaled to the "
                     "atomic constants, calibrated here in the all-electron X2C basis) on the same TDA states")
a = ap.parse_args()

OUT = os.path.join(ROOT, "paper", "references", "results", f"{a.name}.soc_benchmark.json")
res = {"name": a.name, "method": "sf-X2C-1e %s, TDA, Breit-Pauli SOMF (no picture change)" % a.xc.upper(),
       "light_basis": a.light_basis, "heavy_basis": a.heavy_basis, "timings_s": {}}


def save():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(res, f, indent=1)


def heavy_basis(el):
    b = a.heavy_basis
    if b.startswith("unc-"):
        return gto.uncontract(gto.basis.load(b[4:], el))
    return gto.basis.load(b, el)


def basis_for(elements):
    return {el: (heavy_basis(el) if el in ("Br", "I") else a.light_basis) for el in elements}


def somf_ao(mol, dm):
    """Three Cartesian components (real antisymmetric, Hartree) of h_SO, and its one-electron part alone;
    H_SO = i * sum_k h_k s_k as in psscreen.soc."""
    h1 = mol.intor("int1e_pnucxp", comp=3)
    vj, vk, vk2 = jk.get_jk(mol, [dm, dm, dm], ["ijkl,kl->ij", "ijkl,jk->il", "ijkl,li->kj"],
                            intor="int2e_p1vxp1", comp=3, aosym="s1")
    pref = ALPHA ** 2 / 4 * 2          # sigma = 2 s
    return pref * (h1 + vj - 1.5 * (vk + vk2)), pref * h1


def atom_check(el):
    spin = {"C": 2, "N": 3, "B": 1, "F": 1, "Br": 1, "I": 1}[el]
    mol = gto.M(atom=f"{el} 0 0 0", basis=basis_for([el]), spin=spin, verbose=0)
    mf = scf.ROHF(mol).x2c()
    mf.kernel()
    dm = mf.make_rdm1()
    dm = dm[0] + dm[1]
    C, e = mf.mo_coeff, mf.mo_energy
    S = mol.intor("int1e_ovlp")
    pmask = np.array([lab[2].endswith("p") for lab in mol.ao_labels(fmt=False)])
    pchar = np.einsum("pi,pq,qi->pi", C, S, C)[pmask].sum(0)
    homo = int(np.ceil(mol.nelectron / 2)) - 1
    cand = sorted([i for i in range(C.shape[1]) if pchar[i] > 0.9], key=lambda i: abs(e[i] - e[homo]))[:3]
    out = {}
    for tag, h in zip(("somf", "bp1e"), somf_ao(mol, dm)):
        Ip = np.einsum("pi,kpq,qj->kij", C[:, cand], h, C[:, cand])
        out[tag] = float(np.sqrt((Ip ** 2).sum() / 6.0) * CM1_PER_HARTREE)   # zeta, same invariant as psscreen.soc
    out["exp"] = ZETA_EXP_CM.get(el)
    return out


t0 = time.time()
geo = os.path.join(ROOT, "validation", "runs", f"{a.name}.xyz")
if not os.path.exists(geo):
    geo = os.path.join(ROOT, "validation", "runs", "set3", f"{a.name}_c0.xyz")
sym, xyz = read_xyz(geo)
sym = [s.capitalize() for s in sym]
elements = sorted(set(sym) - {"H"})
res["atom_zeta_cm"] = {el: atom_check(el) for el in elements}
print("atomic zeta (cm-1):", json.dumps(res["atom_zeta_cm"]), flush=True)
res["timings_s"]["atoms"] = time.time() - t0
save()
if a.atom_check_only:
    sys.exit(0)

mol = gto.M(atom=to_pyscf_atom(sym, np.asarray(xyz, float)), basis=basis_for(sorted(set(sym))), verbose=4,
            max_memory=int(os.environ.get("PSS_MAX_MEMORY_MB", 40000)))
res["nao"] = mol.nao
mf = dft.RKS(mol, xc=a.xc).x2c().density_fit()
mf.conv_tol = 1e-9
mf.kernel()
assert mf.converged, "SCF did not converge"
res["timings_s"]["scf"] = time.time() - t0
nocc = mol.nelectron // 2


def tda(singlet, n):
    td = tdscf.TDA(mf)
    td.singlet = singlet
    td.nstates = n
    td.conv_tol = 1e-6
    td.kernel()
    X = []
    for x, _ in td.xy:
        x = np.asarray(x).reshape(nocc, -1)
        X.append(x / np.linalg.norm(x))
    return np.asarray(td.e) * EV, X


eS, xS = tda(True, a.nsing)
eT, xT = tda(False, a.ntrip)
res["E_singlets_eV"], res["E_triplets_eV"] = eS.tolist(), eT.tolist()
res["timings_s"]["tda"] = time.time() - t0
save()

def effective_ao():
    """sum_A s_A * (p x (1/r_A) p) with s_A from the X2C free atom in this basis (same invariant as psscreen.soc)."""
    from psscreen.soc import _raw_atom_integrals
    h = np.zeros((3, mol.nao, mol.nao))
    scales = {}
    for ia in range(mol.natm):
        el = mol.atom_pure_symbol(ia)
        if el not in ZETA_EXP_CM:
            continue
        if el not in scales:
            spin = {"C": 2, "N": 3, "B": 1, "F": 1, "Br": 1, "I": 1, "Cl": 1}[el]
            am = gto.M(atom=f"{el} 0 0 0", basis=basis_for([el]), spin=spin, verbose=0)
            amf = scf.ROHF(am).x2c()
            amf.kernel()
            C, e = amf.mo_coeff, amf.mo_energy
            S = am.intor("int1e_ovlp")
            pmask = np.array([lab[2].endswith("p") for lab in am.ao_labels(fmt=False)])
            pchar = np.einsum("pi,pq,qi->pi", C, S, C)[pmask].sum(0)
            homo = int(np.ceil(am.nelectron / 2)) - 1
            cand = sorted([i for i in range(C.shape[1]) if pchar[i] > 0.9], key=lambda i: abs(e[i] - e[homo]))[:3]
            Ip = np.einsum("pi,kpq,qj->kij", C[:, cand], _raw_atom_integrals(am, 0), C[:, cand])
            scales[el] = (ZETA_EXP_CM[el] / CM1_PER_HARTREE) / np.sqrt((Ip ** 2).sum() / 6.0)
        h += scales[el] * _raw_atom_integrals(mol, ia)
    res["effective_scales"] = {k: float(v) for k, v in scales.items()}
    return h


h_somf, h_bp1e = somf_ao(mol, mf.make_rdm1())
res["timings_s"]["somf"] = time.time() - t0
ops = [("somf", h_somf), ("bp1e", h_bp1e)]
if a.effective:
    ops.append(("eff_ae", effective_ao()))
for tag, h in ops:
    blocks = _spin_blocks(1j * np.einsum("pi,kpq,qj->kij", mf.mo_coeff, h, mf.mo_coeff))
    V = np.zeros((len(xS), len(xT)))
    for m, Xs in enumerate(xS):
        for n, Xt in enumerate(xT):
            V[m, n] = np.sqrt(sum(abs(_me(_state("S", Xs), _state(k, Xt), blocks, nocc)) ** 2
                                  for k in ("T0", "T+", "T-"))) * CM1_PER_HARTREE
    res[f"V_S_T_cm_{tag}"] = V.tolist()
    print(tag, "V(S1,Tn) cm-1:", np.round(V[0], 3), flush=True)

pipe = os.path.join(os.path.dirname(geo), os.path.basename(geo).replace(".xyz", ".json"))
if os.path.exists(pipe):
    p = json.load(open(pipe))
    res["pipeline"] = {"E_S1_eV": p.get("E_S1_eV"), "E_triplets_eV": p.get("E_triplets_eV"),
                       "SOC_S1_Tn_cm": p.get("SOC_S1_Tn_cm"), "raw_tddft": p.get("raw_tddft")}
res["timings_s"]["total"] = time.time() - t0
save()
print(a.name, "done", round(time.time() - t0), "s")
