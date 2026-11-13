"""Corrección ΔSCF de S1 y T1 (estados HOMO->LUMO) a coste de estado fundamental.

Por qué: TDDFT sobreestima mucho (0.4-0.8 eV) el S1 de BODIPY y otros colorantes tipo cianina
(el "problema de las cianinas"), mientras que T_n sale razonable. Eso desplaza S1 artificialmente
hacia los tripletes superiores y crea canales ISC espurios. Los métodos ΔSCF/ROKS no sufren ese
artefacto en estas moléculas.

- T1: UKS con M_S = 1 (fundamental del triplete) en la geometría S0 -> E_T (vertical).
- S1: determinante mixto (HOMO_alfa -> LUMO_alfa, M_S = 0) por MOM, y la fórmula de suma
  de Ziegler para quitar la contaminación de espín: E_S = 2 E_mix - E_T.
"""
from __future__ import annotations
import numpy as np
from pyscf import dft, scf

HARTREE_EV = 27.211386245988


def _uks(mol, xc, grid_level, dm0=None, occ=None, mo0=None):
    mf = dft.UKS(mol, xc=xc).density_fit(); mf.grids.level = grid_level
    mf.conv_tol = 1e-8; mf.max_cycle = 150
    if occ is not None:
        mf = scf.addons.mom_occ(mf, mo0, occ)
    mf.kernel(dm0)
    return mf


def delta_scf(mol, xc="b3lyp", grid_level=2, e0=None, mf0=None):
    """Devuelve dict con energías verticales ΔSCF (eV) de S1 y T1 y diagnósticos."""
    if mf0 is None:
        mf0 = dft.RKS(mol, xc=xc).density_fit(); mf0.grids.level = grid_level; mf0.kernel()
    e0 = mf0.e_tot
    nocc = mol.nelectron // 2
    mo = (mf0.mo_coeff, mf0.mo_coeff)
    # triplete M_S=1: alfa ocupa hasta LUMO, beta hasta HOMO-1
    occ_a = np.zeros(len(mf0.mo_energy)); occ_a[:nocc + 1] = 1
    occ_b = np.zeros(len(mf0.mo_energy)); occ_b[:nocc - 1] = 1
    molT = mol.copy(); molT.spin = 2; molT.build(False, False)
    dmT = (mf0.mo_coeff * occ_a) @ mf0.mo_coeff.T, (mf0.mo_coeff * occ_b) @ mf0.mo_coeff.T
    mfT = _uks(molT, xc, grid_level, dm0=np.array(dmT), occ=(occ_a, occ_b), mo0=mo)
    # determinante mixto M_S=0: alfa HOMO->LUMO, beta en el fundamental
    occ_a2 = np.zeros(len(mf0.mo_energy)); occ_a2[:nocc - 1] = 1; occ_a2[nocc] = 1
    occ_b2 = np.zeros(len(mf0.mo_energy)); occ_b2[:nocc] = 1
    dmM = (mf0.mo_coeff * occ_a2) @ mf0.mo_coeff.T, (mf0.mo_coeff * occ_b2) @ mf0.mo_coeff.T
    mfM = _uks(mol, xc, grid_level, dm0=np.array(dmM), occ=(occ_a2, occ_b2), mo0=mo)
    ET = (mfT.e_tot - e0) * HARTREE_EV
    EM = (mfM.e_tot - e0) * HARTREE_EV
    ES = 2 * EM - ET
    return {"E_S1_dscf": ES, "E_T1_dscf": ET, "E_mixed": EM,
            "S2_T": float(mfT.spin_square()[0]), "S2_mixed": float(mfM.spin_square()[0]),
            "converged": bool(mfT.converged and mfM.converged)}
