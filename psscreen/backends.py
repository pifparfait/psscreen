"""Backends de estados excitados.

Contrato común (ExcitedStateResult): energías verticales de singletes y tripletes (eV),
fuerzas de oscilador, y amplitudes TDA X[nocc,nvir] en la base de MOs, más los MOs.
Con eso, el módulo SOC y el de tasas funcionan igual sea cual sea el backend.

- PySCFTDA: TDDFT (RPA) para singletes y TDA para tripletes, en ventana de orbitales (tdwin.py).
- ActiveSpaceBackend: hueco para un backend cuántico (QPE/QSP sobre un espacio activo).
  Exporta el Hamiltoniano activo (FCIDUMP) y, como referencia clásica del mismo problema,
  lo resuelve exactamente con CASCI (lo que un ordenador cuántico tolerante a fallos haría
  con QPE) y opcionalmente añade correlación dinámica con NEVPT2.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from abc import ABC, abstractmethod
import numpy as np

HARTREE_EV = 27.211386245988


@dataclass
class ExcitedStateResult:
    mol: object
    mo_coeff: np.ndarray
    nocc: int
    e_singlets: np.ndarray            # eV
    f_osc: np.ndarray
    e_triplets: np.ndarray            # eV
    x_singlets: list = field(default_factory=list)   # cada uno [nocc, nvir], norma 1
    x_triplets: list = field(default_factory=list)
    meta: dict = field(default_factory=dict)


class ExcitedStateBackend(ABC):
    name = "abstract"

    @abstractmethod
    def run(self, atom: str, charge: int = 0) -> ExcitedStateResult: ...


_HEAVY = {"Rb", "Sr", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Pt", "Au", "Hg", "Pb", "Bi", "Ir", "Ru", "Pd", "Os"}


def _ecp(atom: str, basis: str):
    """ECP solo para elementos pesados (def2: Z>36)."""
    els = {l.split()[0] for l in atom.strip().splitlines()}
    return {e: basis for e in els if e in _HEAVY} or None


class PySCFTDA(ExcitedStateBackend):
    name = "pyscf-tda"

    def __init__(self, xc="b3lyp", basis="def2-svp", nstates=6, solvent=None,
                 density_fit=True, verbose=0, max_memory=12000, ewin=10.0, rpa=True,
                 max_occ=60, max_vir=80, nsinglets=2, ewin_singlets=60.0, grid_level=1,
                 scf_grid_level=2, delta_scf=True):
        self.xc, self.basis, self.nstates = xc, basis, nstates
        self.ewin, self.rpa, self.max_occ, self.max_vir = ewin, rpa, max_occ, max_vir
        self.nsinglets, self.ewin_singlets, self.grid_level = nsinglets, ewin_singlets, grid_level
        self.scf_grid_level = scf_grid_level
        self.delta_scf = delta_scf
        self.solvent, self.density_fit = solvent, density_fit
        self.verbose, self.max_memory = verbose, max_memory

    def run(self, atom, charge=0):
        from pyscf import gto, dft, tdscf
        mol = gto.M(atom=atom, basis=self.basis, ecp=_ecp(atom, self.basis), charge=charge,
                    verbose=self.verbose, max_memory=self.max_memory)
        mf = dft.RKS(mol, xc=self.xc)
        mf.grids.level = self.scf_grid_level
        if self.density_fit:
            mf = mf.density_fit()
        if self.solvent:
            raise NotImplementedError("PCM aún no conectado al TDDFT en ventana")
        mf.conv_tol = 1e-9
        mf.kernel()
        if not mf.converged:
            raise RuntimeError("SCF no convergió")
        from .tdwin import window_td, rpa_singlets_davidson
        import time
        nocc = mol.nelectron // 2
        t0 = time.time()
        # tripletes (y singletes TDA de control) en ventana: convergen rápido con la ventana
        r = window_td(mf, ewin_ev=self.ewin, nstates=self.nstates, rpa_singlets=False,
                      max_occ=self.max_occ, max_vir=self.max_vir, grid_level=self.grid_level)
        t1 = time.time()
        # singletes RPA en espacio (casi) completo: el estado brillante necesita la respuesta
        # de los virtuales altos (en ventanas pequeñas sale ~0.2-0.4 eV demasiado alto)
        s = rpa_singlets_davidson(mf, nstates=self.nsinglets, ewin_ev=self.ewin_singlets,
                                  grid_level=0 if self.grid_level is None else min(self.grid_level, 1))
        if not s["converged"]:
            raise RuntimeError("Davidson RPA no convergió")
        t2 = time.time()
        dscf = None
        if self.delta_scf:
            from .deltascf import delta_scf
            dscf = delta_scf(mol, xc=self.xc, grid_level=self.scf_grid_level, mf0=mf)
        t3 = time.time()
        return ExcitedStateResult(mol=mol, mo_coeff=mf.mo_coeff, nocc=nocc,
                                  e_singlets=s["e_S"] * HARTREE_EV, f_osc=s["f"],
                                  e_triplets=r["e_T"] * HARTREE_EV,
                                  x_singlets=s["x_S"], x_triplets=r["x_T"],
                                  meta={"backend": self.name, "xc": self.xc, "basis": self.basis,
                                        "e_scf": mf.e_tot, "nao": mol.nao, "ewin_eV": self.ewin,
                                        "n_occ_win": len(r["occ"]), "n_vir_win": len(r["vir"]),
                                        "singlets": f"RPA Davidson, ventana {self.ewin_singlets} eV, npair {s['npair']}",
                                        "triplets": f"TDA explícito, ventana {self.ewin} eV",
                                        "S_TDA_window_eV": (r["e_S"] * HARTREE_EV).tolist(),
                                        "t_triplets_s": t1 - t0, "t_singlets_s": t2 - t1,
                                        "delta_scf": dscf, "t_dscf_s": t3 - t2})


class ActiveSpaceBackend(ExcitedStateBackend):
    """Hueco para el backend cuántico.

    `solver` decide quién diagonaliza el Hamiltoniano del espacio activo:
      - "casci" (por defecto): diagonalización exacta clásica = resultado ideal de QPE.
      - un callable externo f(h1, h2, ncas, nelecas, nroots, spin) -> energías (Ha),
        p.ej. un estimador QPE/QSP en hardware o simulador.
    `nevpt2=True` añade correlación dinámica (SC-NEVPT2) sobre CASCI, que es lo que el
    artículo de Xanadu no incluye.
    Devuelve solo energías (sin amplitudes TDA), así que el SOC/ISC se toma del backend TDA
    salvo que se proporcione un modelo SOC propio.
    """
    name = "active-space"

    def __init__(self, ncas=6, nelecas=6, basis="def2-svp", nroots=3, solver="casci",
                 nevpt2=False, fcidump_path=None, verbose=0):
        self.ncas, self.nelecas, self.basis, self.nroots = ncas, nelecas, basis, nroots
        self.solver, self.nevpt2, self.fcidump_path = solver, nevpt2, fcidump_path
        self.verbose = verbose

    def run(self, atom, charge=0):
        from pyscf import gto, scf, mcscf, fci
        mol = gto.M(atom=atom, basis=self.basis, ecp=_ecp(atom, self.basis), charge=charge, verbose=self.verbose)
        mf = scf.RHF(mol).density_fit().run()
        res = {}
        for mult, spin in (("S", 0), ("T", 2)):
            mc = mcscf.CASCI(mf, self.ncas, self.nelecas)
            mc.fcisolver = fci.direct_spin1.FCI(mol)
            mc.fcisolver.nroots = self.nroots + (1 if spin == 0 else 0)
            fci.addons.fix_spin_(mc.fcisolver, ss=spin / 2 * (spin / 2 + 1))
            if self.fcidump_path and spin == 0:
                from pyscf.tools import fcidump
                h1, ecore = mc.get_h1eff(); h2 = mc.get_h2eff()
                fcidump.from_integrals(self.fcidump_path, h1, h2, self.ncas, self.nelecas, ecore)
            if callable(self.solver):
                h1, ecore = mc.get_h1eff(); h2 = mc.get_h2eff()
                e = np.asarray(self.solver(h1, h2, self.ncas, self.nelecas, mc.fcisolver.nroots, spin)) + ecore
            else:
                mc.kernel(); e = np.atleast_1d(mc.e_tot)
                if self.nevpt2:
                    from pyscf import mrpt
                    e = np.array([mc.e_tot[i] + mrpt.NEVPT(mc, root=i).kernel() for i in range(len(e))])
            res[mult] = np.sort(e)
        e0 = res["S"][0]
        return ExcitedStateResult(mol=mol, mo_coeff=mf.mo_coeff, nocc=mol.nelectron // 2,
                                  e_singlets=(res["S"][1:] - e0) * HARTREE_EV,
                                  f_osc=np.full(len(res["S"]) - 1, np.nan),
                                  e_triplets=(res["T"] - e0) * HARTREE_EV,
                                  meta={"backend": self.name, "ncas": self.ncas, "nelecas": self.nelecas,
                                        "nevpt2": self.nevpt2})
