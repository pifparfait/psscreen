"""Controlled pi active space for the BODIPY dyes: dipyrrin 2p_z AVAS selection, SA-CASSCF, CASCI and SC-NEVPT2.

The earlier multireference runs (run_sa_nevpt2.py, run_cas_sweep.py) take the active orbitals from the canonical RHF
orbitals around the Fermi level, so the same (ne,no) label can hold different kinds of orbitals in different dyes
(pi, iodine lone pairs, sigma). This script fixes the active space by a reproducible rule:

  1. Find the dipyrrin pi core: the two pyrrole rings bonded to boron through N and the meso atom (C, or N in
     aza-BODIPYs) that bridges them. That is 11 atoms with 12 pi electrons, the (12e,11o) space of Zhou et al.
  2. Rotate the molecule so that the best plane through those atoms is the xy plane.
  3. AVAS-type selection: project the occupied and the virtual RHF orbitals separately onto the MINAO 2p_z functions
     of the 11 core atoms and keep the 6 occupied and 5 virtual orbitals with the largest projection.
     The projection eigenvalues of the kept and of the first rejected orbitals are stored.
  4. CASCI singlets/triplets and SC-NEVPT2 on these orbitals (no optimisation): compare with the RHF-window runs.
  5. SA-CASSCF over the lowest singlets (equal weights), then CASCI singlets/triplets and SC-NEVPT2.
It also stores alpha_SO (paper/quantum/alpha_so.py) in both spaces. For every set of orbitals the script stores the composition of each active orbital (weight on the core 2p_z space,
Mulliken weight on iodine/bromine, on the core atoms and on the rest) and the natural occupation numbers of the
state-averaged and of each root's density.

Usage (repository root): python paper/references/run_pi_cas.py NAME [--nsing 3] [--ntrip 3] [--no-nevpt2]
                         [--basis def2-tzvp] [--sa-mix] [--skip-avas-block] [--select-only]
Environment: OMP_NUM_THREADS, PSS_MAX_MEMORY_MB (default 40000).
Writes paper/references/results/NAME.pi_cas_12_11[_<basis>][_samix].json, updated after every step.
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
NOCC_ACT, NVIR_ACT = 6, 5          # (12e, 11o)

ap = argparse.ArgumentParser()
ap.add_argument("name")
ap.add_argument("--nsing", type=int, default=3)
ap.add_argument("--ntrip", type=int, default=3)
ap.add_argument("--max-macro", type=int, default=100)
ap.add_argument("--no-nevpt2", action="store_true")
ap.add_argument("--basis", default="def2-svp")
ap.add_argument("--sa-mix", action="store_true",
                help="state-average over the singlets and the two lowest triplets (equal weights) instead of singlets only")
ap.add_argument("--skip-avas-block", action="store_true", help="skip CASCI/NEVPT2 on the unoptimized AVAS orbitals")
ap.add_argument("--select-only", action="store_true", help="stop after the orbital selection (quick check)")
a = ap.parse_args()

geo = os.path.join(ROOT, "validation", "runs", f"{a.name}.xyz")
if not os.path.exists(geo):
    geo = os.path.join(ROOT, "validation", "runs", "set3", f"{a.name}_c0.xyz")
sym, xyz = read_xyz(geo)
sym = [s.capitalize() for s in sym]
xyz = np.asarray(xyz, float)


# ---------------------------------------------------------------- 1) dipyrrin core
def bonds(sym, xyz):
    rcov = {"H": 0.32, "B": 0.85, "C": 0.75, "N": 0.71, "O": 0.63, "F": 0.64, "S": 1.03, "Cl": 1.0, "Br": 1.14,
            "I": 1.33}
    n = len(sym)
    d = np.linalg.norm(xyz[:, None] - xyz[None], axis=-1)
    adj = [set() for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if d[i, j] < 1.25 * (rcov.get(sym[i], 1.0) + rcov.get(sym[j], 1.0)):
                adj[i].add(j)
                adj[j].add(i)
    return adj


def five_ring(adj, start, allowed):
    """The 5-membered ring through `start` whose atoms are all in `allowed`."""
    def dfs(path):
        if len(path) == 5:
            return path if start in adj[path[-1]] else None
        for nb in adj[path[-1]]:
            if nb in allowed and nb not in path:
                r = dfs(path + [nb])
                if r:
                    return r
        return None
    return dfs([start])


adj = bonds(sym, xyz)
boron = [i for i, s in enumerate(sym) if s == "B"]
assert len(boron) == 1, "expected one boron"
cn = {i for i, s in enumerate(sym) if s in ("C", "N")}
n_b = [j for j in adj[boron[0]] if sym[j] == "N"]
assert len(n_b) == 2, "expected two N bonded to B"
rings = [five_ring(adj, n, cn) for n in n_b]
assert all(rings), "pyrrole ring not found"
ring_atoms = set(rings[0]) | set(rings[1])
meso = [k for k in cn - ring_atoms
        if any(x in rings[0] for x in adj[k]) and any(x in rings[1] for x in adj[k])]
assert len(meso) == 1, f"meso atom not unique: {meso}"
core = sorted(ring_atoms | set(meso))
assert len(core) == 11

# ---------------------------------------------------------------- 2) rotate core plane to xy
c0 = xyz[core].mean(0)
_, _, vt = np.linalg.svd(xyz[core] - c0)
R = vt                                            # rows: in-plane x, in-plane y, normal
xyz_r = (xyz - c0) @ R.T
planarity = float(np.abs(xyz_r[core, 2]).max())
atom = to_pyscf_atom(sym, xyz_r)
mol = gto.M(atom=atom, basis=a.basis, ecp=_ecp(atom, a.basis), verbose=4, unit="Angstrom",
            max_memory=float(os.environ.get("PSS_MAX_MEMORY_MB", 40000)))

outdir = os.path.join(ROOT, "paper", "references", "results")
os.makedirs(outdir, exist_ok=True)
tag = ("" if a.basis == "def2-svp" else "_" + a.basis.replace("-", "")) + ("_samix" if a.sa_mix else "")
fout = os.path.join(outdir, f"{a.name}.pi_cas_12_11{tag}.json")
out = {"name": a.name, "geometry": os.path.relpath(geo, ROOT), "basis": a.basis, "cas": [12, 11],
       "sa_weights": "singlets + 2 triplets, equal" if a.sa_mix else "singlets only, equal",
       "selection": "AVAS-type projection onto MINAO 2p_z of the 11 dipyrrin atoms, 6 occupied + 5 virtual",
       "core_atoms": [f"{sym[i]}{i}" for i in core], "core_max_out_of_plane_A": planarity,
       "nsing": a.nsing, "ntrip": a.ntrip, "timings_s": {}}


def save():
    json.dump(out, open(fout, "w"), indent=1)


t0 = time.time()
mf = scf.RHF(mol).density_fit().run()
out["timings_s"]["rhf"] = time.time() - t0

# ---------------------------------------------------------------- 3) AVAS-type selection
pmol = gto.M(atom=to_pyscf_atom([sym[i] for i in core], xyz_r[core]), basis="minao", unit="Angstrom", verbose=0)
pz = [k for k, lab in enumerate(pmol.ao_labels()) if lab.split()[2] == "2pz"]
assert len(pz) == 11
S = mf.get_ovlp()
S12 = gto.intor_cross("int1e_ovlp", mol, pmol)[:, pz]
S22 = pmol.intor("int1e_ovlp")[np.ix_(pz, pz)]
PS = S12 @ np.linalg.solve(S22, S12.T)            # AO representation of the projector onto core 2p_z (with S)
nocc = mol.nelectron // 2
C = mf.mo_coeff
Co, Cv = C[:, :nocc], C[:, nocc:]
wo, uo = np.linalg.eigh(Co.T @ PS @ Co)
wv, uv = np.linalg.eigh(Cv.T @ PS @ Cv)
io, iv = np.argsort(-wo), np.argsort(-wv)
Co_r, Cv_r = Co @ uo[:, io], Cv @ uv[:, iv]
out["avas_occ_eigs"] = [float(x) for x in wo[io][:NOCC_ACT + 4]]
out["avas_vir_eigs"] = [float(x) for x in wv[iv][:NVIR_ACT + 4]]
print(a.name, "AVAS occ eigs", np.round(wo[io][:NOCC_ACT + 3], 3), "vir eigs", np.round(wv[iv][:NVIR_ACT + 3], 3),
      flush=True)


def semicanon(Cblk):
    f = Cblk.T @ mf.get_fock() @ Cblk
    e, u = np.linalg.eigh(f)
    return Cblk @ u


core_mo = semicanon(Co_r[:, NOCC_ACT:])
act_mo = semicanon(np.hstack([Co_r[:, :NOCC_ACT], Cv_r[:, :NVIR_ACT]]))
vir_mo = semicanon(Cv_r[:, NVIR_ACT:])
mo_avas = np.hstack([core_mo, act_mo, vir_mo])
ncore = core_mo.shape[1]

heavy = [k for k, lab in enumerate(mol.ao_labels()) if lab.split()[1] in ("I", "Br")]
core_ao = [k for k, lab in enumerate(mol.ao_labels()) if int(lab.split()[0]) in core]


def composition(mo, first=None):
    """Per active orbital: core 2p_z weight, and Mulliken weights on core atoms, on I/Br, and elsewhere."""
    first = ncore if first is None else first
    cas = mo[:, first:first + 11]
    res = []
    for k in range(11):
        c = cas[:, k]
        mull = c * (S @ c)
        w_core = float(mull[core_ao].sum())
        w_heavy = float(mull[heavy].sum()) if heavy else 0.0
        res.append({"pz_core": float(c @ PS @ c), "mulliken_core_atoms": w_core, "mulliken_I_Br": w_heavy,
                    "mulliken_rest": 1.0 - w_core - w_heavy})
    return res


out["composition_avas"] = composition(mo_avas)
# the canonical RHF window (HOMO-5 ... LUMO+4) used by run_sa_nevpt2.py and run_cas_sweep.py, for comparison
out["composition_canonical_window"] = composition(mf.mo_coeff, nocc - NOCC_ACT)
out["sum_pz_core"] = {"avas": sum(c["pz_core"] for c in out["composition_avas"]),
                      "canonical_window": sum(c["pz_core"] for c in out["composition_canonical_window"])}


def alpha_so_cm(mo, first):
    """Block-encoding norm of the one-body SOC operator in the 11 orbitals starting at `first` (alpha_so.py)."""
    from psscreen.soc import soc_mo_integrals
    sys.path.insert(0, os.path.join(ROOT, "paper", "quantum"))
    from alpha_so import HA_CM, spin_orbital_matrix
    h, _ = soc_mo_integrals(mol, mo, a.basis)
    idx = list(range(first, first + 11))
    M = spin_orbital_matrix(h[:, idx][:, :, idx])
    return 0.5 * float(np.abs(np.linalg.eigvalsh(0.5 * (M + M.conj().T))).sum()) * HA_CM


try:
    out["alpha_so_cm_hf_orbitals"] = {"avas": alpha_so_cm(mo_avas, ncore),
                                      "canonical_window": alpha_so_cm(mf.mo_coeff, nocc - NOCC_ACT)}
    print(a.name, "alpha_SO (cm-1)", out["alpha_so_cm_hf_orbitals"], flush=True)
except Exception as exc:          # never lose the CAS runs to this side calculation
    out["alpha_so_error"] = repr(exc)
save()
if a.select_only:
    print(a.name, "selection only; core", out["core_atoms"], "planarity", round(planarity, 3), flush=True)
    sys.exit(0)


def casci(mo, spin, nroots):
    c = mcscf.CASCI(mf, 11, 12)
    if spin == 0:
        c.fcisolver = fci.direct_spin0.FCI(mol)
    else:
        c.fcisolver = fci.direct_spin1.FCI(mol)
        c.nelecas = ((12 + spin) // 2, (12 - spin) // 2)
        c.fix_spin_(ss=spin / 2 * (spin / 2 + 1))
    c.fcisolver.nroots = nroots
    c.kernel(mo)
    return c


def natocc(c, nroots, spin):
    """Natural occupations of each root, and of the equal-weight state-averaged density."""
    occ, dms = [], []
    for i in range(nroots):
        ci = c.ci[i] if nroots > 1 else c.ci
        dm = c.fcisolver.make_rdm1(ci, 11, c.nelecas)
        dms.append(dm)
        occ.append([float(x) for x in np.sort(np.linalg.eigvalsh(dm))[::-1]])
    sa = [float(x) for x in np.sort(np.linalg.eigvalsh(np.mean(dms, axis=0)))[::-1]]
    return occ, sa


def block(tag, mo):
    t = time.time()
    cs = casci(mo, 0, a.nsing)
    ct = casci(mo, 2, a.ntrip)
    E0 = cs.e_tot[0]
    r = {"E_CASCI_singlets_eV": [(e - E0) * EV for e in cs.e_tot],
         "E_CASCI_triplets_eV": [(e - E0) * EV for e in ct.e_tot],
         }
    r["natocc_singlets"], r["natocc_sa_singlets"] = natocc(cs, a.nsing, 0)   # SA: occupations of the averaged density
    r["natocc_triplets"], r["natocc_sa_triplets"] = natocc(ct, a.ntrip, 2)
    out[tag] = r
    out["timings_s"][tag + "_casci"] = time.time() - t
    print(a.name, tag, "CASCI S:", np.round(r["E_CASCI_singlets_eV"], 3), "T:", np.round(r["E_CASCI_triplets_eV"], 3),
          flush=True)
    print(a.name, tag, "SA natocc:", np.round(r["natocc_sa_singlets"], 3), flush=True)
    save()
    if a.no_nevpt2:
        return
    t = time.time()
    try:
        es = [cs.e_tot[i] + mrpt.NEVPT(cs, root=i).kernel() for i in range(a.nsing)]
        r["E_NEVPT2_singlets_eV"] = [(e - es[0]) * EV for e in es]
        save()
        et = [ct.e_tot[i] + mrpt.NEVPT(ct, root=i).kernel() for i in range(a.ntrip)]
        r["E_NEVPT2_triplets_eV"] = [(e - es[0]) * EV for e in et]
        print(a.name, tag, "NEVPT2 S:", np.round(r["E_NEVPT2_singlets_eV"], 3), "T:",
              np.round(r["E_NEVPT2_triplets_eV"], 3), flush=True)
    except Exception as exc:
        r["nevpt2_error"] = repr(exc)
        print(a.name, tag, "NEVPT2 failed:", exc, flush=True)
    out["timings_s"][tag + "_nevpt2"] = time.time() - t
    save()


# ---------------------------------------------------------------- 4) on the AVAS orbitals
if not a.skip_avas_block:
    block("avas_orbitals", mo_avas)

# ---------------------------------------------------------------- 5) SA-CASSCF from the AVAS orbitals
t = time.time()
mc = mcscf.CASSCF(mf, 11, 12).density_fit()
if a.sa_mix:
    fs = fci.direct_spin1.FCI(mol)
    fs.spin, fs.nroots = 0, a.nsing
    fs = fci.addons.fix_spin_(fs, ss=0)
    ft = fci.direct_spin1.FCI(mol)
    ft.spin, ft.nroots = 2, 2
    ft = fci.addons.fix_spin_(ft, ss=2)
    nst = a.nsing + 2
    mc = mcscf.addons.state_average_mix_(mc, [fs, ft], [1.0 / nst] * nst)
else:
    mc.fcisolver = fci.direct_spin0.FCI(mol)
    mc.fcisolver.nroots = a.nsing
    mc = mcscf.state_average_(mc, [1.0 / a.nsing] * a.nsing)
mc.max_cycle_macro = a.max_macro
mc.kernel(mo_avas)
out["sa_casscf_converged"] = bool(mc.converged)
out["timings_s"]["sa_casscf"] = time.time() - t
out["composition_sa_casscf"] = composition(mc.mo_coeff)
save()
block("sa_casscf_orbitals", mc.mo_coeff)
out["timings_s"]["total"] = time.time() - t0
save()
print(a.name, "done", round(time.time() - t0), "s", flush=True)
