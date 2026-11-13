"""Effective one-electron spin-orbit coupling (SOC), calibrated against atomic spectroscopy.

Operator: H_SO = sum_A s_A * L_A . s , with L_A ~ (r_A x p)/r_A^3 (PySCF int1e_prinvxp integral).
The per-element factor s_A is NOT taken from Z_eff tables: it is calibrated so that the isolated
atom, in the SAME basis/ECP as the molecule, reproduces the experimental zeta_np constant of its
valence shell. The constant thus absorbs two-electron screening, scalar relativistic effects and
the shape of the ECP pseudo-orbitals (a Koseki-type approach, but self-calibrated).

Matrix elements <S_m|H_SO|T_n,M> between spin-adapted TDA states with a closed-shell reference,
summed over M = -1,0,+1 (what enters the golden rule).
"""
from __future__ import annotations
import numpy as np
from functools import lru_cache

CM1_PER_HARTREE = 219474.6313632

# Experimental zeta_np (cm^-1) of the valence shell of the neutral atom.
# Halogens: 2/3 of the 2P3/2-2P1/2 splitting (NIST ASD): F 404.1, Cl 882.4, Br 3685.2, I 7603.2 cm^-1.
# B: 2/3 * 15.29 (2P). C, N, O: standard table values (Montalti et al., Handbook of
# Photochemistry, 3rd ed., table of atomic spin-orbit constants).
ZETA_EXP_CM = {"B": 10.2, "C": 28.0, "N": 76.0, "O": 151.0, "F": 269.4, "S": 382.0,
               "Cl": 588.3, "Br": 2456.8, "I": 5068.8}
VALENCE_N = {"B": 2, "C": 2, "N": 2, "O": 2, "F": 2, "S": 3, "Cl": 3, "Br": 4, "I": 5}
SPIN_ATOM = {"B": 1, "C": 2, "N": 3, "O": 2, "F": 1, "S": 2, "Cl": 1, "Br": 1, "I": 1}


def _raw_atom_integrals(mol, ia):
    with mol.with_rinv_origin(mol.atom_coord(ia)):
        return mol.intor("int1e_prinvxp", comp=3)   # real, antisymmetric, atomic units


@lru_cache(maxsize=None)
def calibrate_element(elem: str, basis: str) -> float:
    """Return s_A (Hartree per unit integral) such that zeta_calc(atom) = zeta_exp."""
    from pyscf import gto, scf
    mol = gto.M(atom=f"{elem} 0 0 0", basis=basis, ecp=({elem: basis} if elem in ("I", "Te", "Xe", "Sb") else None), spin=SPIN_ATOM[elem], verbose=0)
    mf = scf.ROHF(mol)
    mf.kernel()
    C = mf.mo_coeff
    I = _raw_atom_integrals(mol, 0)
    # identify the 3 valence np MOs: highest p character in the valence shell
    labels = mol.ao_labels(fmt=False)
    S = mol.intor("int1e_ovlp")
    pmask = np.array([lab[2].endswith("p") for lab in labels])
    pop = np.einsum("pi,pq,qi->pi", C, S, C)  # Mulliken per AO
    pchar = pop[pmask].sum(0)
    nocc_tot = int(np.ceil(mol.nelectron / 2))
    cand = [i for i in range(C.shape[1]) if pchar[i] > 0.9]
    # take the p triplet with energy closest to the HOMO (valence shell)
    homo_e = mf.mo_energy[nocc_tot - 1]
    cand = sorted(cand, key=lambda i: abs(mf.mo_energy[i] - homo_e))[:3]
    Cp = C[:, cand]
    Ip = np.einsum("pi,kpq,qj->kij", Cp, I, Cp)
    # rotational invariant: sum_k,ij |zeta*l_k|^2 = 6 zeta^2 for l=1
    zeta_unit = np.sqrt((Ip ** 2).sum() / 6.0)
    return (ZETA_EXP_CM[elem] / CM1_PER_HARTREE) / zeta_unit


def soc_mo_integrals(mol, mo_coeff, basis: str) -> np.ndarray:
    """h[k,p,q] (Hartree, complex) in the MO basis; H_SO = sum_k h_k s_k."""
    nao = mol.nao
    h = np.zeros((3, nao, nao))
    scales = {}
    for ia in range(mol.natm):
        el = mol.atom_pure_symbol(ia)
        if el not in ZETA_EXP_CM:
            continue  # H: despreciable
        if el not in scales:
            scales[el] = calibrate_element(el, basis)
        h += scales[el] * _raw_atom_integrals(mol, ia)
    h_mo = np.einsum("pi,kpq,qj->kij", mo_coeff, h, mo_coeff)
    return 1j * h_mo, scales


_PAULI = {"x": np.array([[0, 1], [1, 0]]) / 2, "y": np.array([[0, -1j], [1j, 0]]) / 2,
          "z": np.array([[1, 0], [0, -1]]) / 2}


def _spin_blocks(h_mo):
    """h^{mu,nu}_pq = sum_k h^k_pq <mu|s_k|nu>, mu,nu in {0:alpha,1:beta}."""
    blocks = np.zeros((2, 2) + h_mo.shape[1:], dtype=complex)
    for k, key in enumerate("xyz"):
        for mu in range(2):
            for nu in range(2):
                blocks[mu, nu] += h_mo[k] * _PAULI[key][mu, nu]
    return blocks


def _state(kind, X):
    r2 = 1 / np.sqrt(2)
    if kind == "S":
        return {(0, 0): r2 * X, (1, 1): r2 * X}
    if kind == "T0":
        return {(0, 0): r2 * X, (1, 1): -r2 * X}
    if kind == "T+":   # alpha virtual created, beta occupied destroyed: a†_aα a_iβ
        return {(1, 0): X}
    if kind == "T-":
        return {(0, 1): X}
    raise ValueError(kind)


def _me(bra, ket, blocks, nocc):
    """<bra|H|ket> for combinations of single excitations with a spin flip.
    Determinant (i s1 -> a s2). Formula: d_ij d_s1t1 h^{s2 t2}_ab - d_ab d_s2t2 h^{t1 s1}_ji."""
    val = 0j
    for (s1, s2), cb in bra.items():
        for (t1, t2), ck in ket.items():
            if s1 == t1:
                hvv = blocks[s2, t2][nocc:, nocc:]
                val += np.einsum("ia,ib,ab->", cb.conj(), ck, hvv)
            if s2 == t2:
                hoo = blocks[t1, s1][:nocc, :nocc]
                val -= np.einsum("ia,ja,ji->", cb.conj(), ck, hoo)
    return val


def soc_singlet_triplet(res, h_mo) -> np.ndarray:
    """Matriz V[m,n] = sqrt(sum_M |<S_m|H|T_n,M>|^2) en cm^-1."""
    blocks = _spin_blocks(h_mo)
    nS, nT = len(res.x_singlets), len(res.x_triplets)
    V = np.zeros((nS, nT))
    for m, Xs in enumerate(res.x_singlets):
        bra = _state("S", Xs)
        for n, Xt in enumerate(res.x_triplets):
            tot = sum(abs(_me(bra, _state(k, Xt), blocks, res.nocc)) ** 2 for k in ("T0", "T+", "T-"))
            V[m, n] = np.sqrt(tot) * CM1_PER_HARTREE
    return V


def soc_ground_triplet(res, h_mo) -> np.ndarray:
    """<S0|H|T_n> (cm^-1), useful for phosphorescence / T1->S0 ISC."""
    blocks = _spin_blocks(h_mo)
    out = []
    for Xt in res.x_triplets:
        tot = 0.0
        for k in ("T0", "T+", "T-"):
            ket = _state(k, Xt); v = 0j
            for (t1, t2), ck in ket.items():
                v += np.einsum("ia,ia->", blocks[t1, t2][:res.nocc, res.nocc:], ck)
            tot += abs(v) ** 2
        out.append(np.sqrt(tot) * CM1_PER_HARTREE)
    return np.array(out)
