"""Pipeline de cribado: geometría -> estados excitados -> SOC -> tasas -> puntuación PDT."""
from __future__ import annotations
import json, os, time
import numpy as np
from .geometry import smiles_to_xyz, xtb_optimize, to_pyscf_atom, write_xyz, read_xyz
from .backends import PySCFTDA, ExcitedStateBackend
from .soc import soc_mo_integrals, soc_singlet_triplet, soc_ground_triplet
from .rates import photophysics, MLJParams

EV_NM = 1239.84198
E_O2_SINGLET = 0.98  # eV, energía de 1Δg del O2: T1 debe superarla para fotosensibilizar tipo II


def nm(e):
    return EV_NM / e


def pdt_score(lam_nm, f, phi_isc, E_T1, window=(650, 850)):
    """Puntuación 0-1 orientativa: absorción en ventana terapéutica x rendimiento triplete x T1>1Δg."""
    lo, hi = window
    w = 1.0 if lo <= lam_nm <= hi else float(np.exp(-((min(abs(lam_nm - lo), abs(lam_nm - hi))) / 50.0) ** 2))
    bright = min(1.0, f / 0.3)
    o2 = 1.0 if E_T1 > E_O2_SINGLET + 0.05 else 0.0
    return w * bright * phi_isc * o2


def screen_molecule(name, smiles=None, xyz_path=None, backend: ExcitedStateBackend | None = None,
                    workdir="runs", k_ic=1.5e8, n_refr=1.0, mlj=MLJParams(), corrections=None,
                    window=(650, 850)):
    os.makedirs(workdir, exist_ok=True)
    t0 = time.time()
    geom = os.path.join(workdir, f"{name}.xyz")
    if xyz_path:
        sym, xyz = read_xyz(xyz_path)
    elif os.path.exists(geom):
        sym, xyz = read_xyz(geom)
    else:
        sym, xyz = smiles_to_xyz(smiles)
        sym, xyz, _ = xtb_optimize(sym, xyz)
        write_xyz(geom, sym, xyz, f"{name} GFN2-xTB {smiles}")
    backend = backend or PySCFTDA()
    res = backend.run(to_pyscf_atom(sym, xyz))
    t1 = time.time()
    h_mo, scales = soc_mo_integrals(res.mol, res.mo_coeff, backend.basis)
    V = soc_singlet_triplet(res, h_mo)
    V0 = soc_ground_triplet(res, h_mo)
    # S1 = primer singlete brillante entre los dos más bajos (en BODIPY es S1)
    iS1 = 0
    corr = dict(corrections or {})
    d = res.meta.get("delta_scf")
    if d and d.get("converged") and not corr:
        # S1/T1 de ΔSCF sustituyen a los de TDDFT; el resto de la escalera se desplaza igual
        corr = {"S": d["E_S1_dscf"] - res.e_singlets[0], "T": d["E_T1_dscf"] - res.e_triplets[0],
                "source": "delta-SCF"}
    ph = photophysics(res.e_singlets[iS1], res.f_osc[iS1], res.e_triplets, V[iS1], k_ic=k_ic,
                      n_refr=n_refr, p=mlj, shift_S=corr.get("S", 0.0), shift_T=corr.get("T", 0.0))
    E_S1 = res.e_singlets[iS1] + corr.get("S", 0.0)
    E_T1 = res.e_triplets[0] + corr.get("T", 0.0)
    ph_raw = photophysics(res.e_singlets[iS1], res.f_osc[iS1], res.e_triplets, V[iS1], k_ic=k_ic,
                          n_refr=n_refr, p=mlj)
    out = {
        "name": name, "corrections": {k: (float(v) if not isinstance(v, str) else v) for k, v in corr.items()},
        "raw_tddft": {"E_S1_eV": float(res.e_singlets[iS1]), "lambda_nm": nm(res.e_singlets[iS1]),
                      "E_T1_eV": float(res.e_triplets[0]), "photophysics": ph_raw}, "smiles": smiles, "natoms": len(sym), "meta": res.meta,
        "E_singlets_eV": res.e_singlets.tolist(), "f_osc": res.f_osc.tolist(),
        "E_triplets_eV": res.e_triplets.tolist(),
        "lambda_max_nm": nm(E_S1), "E_S1_eV": float(E_S1), "E_T1_eV": float(E_T1),
        "dE_ST_eV": float(E_S1 - E_T1),
        "SOC_S1_Tn_cm": V[iS1].tolist(), "SOC_S0_Tn_cm": V0.tolist(),
        "soc_scales": scales, "photophysics": ph,
        "proxy_xanadu_sum_V2": float((V[iS1] ** 2).sum()),
        "pdt_score": pdt_score(nm(E_S1), res.f_osc[iS1], ph["phi_ISC"], E_T1, window),
        "wall_s": {"excited": t1 - t0, "total": time.time() - t0},
    }
    with open(os.path.join(workdir, f"{name}.json"), "w") as f:
        json.dump(out, f, indent=1, default=float)
    return out
