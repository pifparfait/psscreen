"""Backend cuántico (PennyLane, simulador) sobre el conjunto de validación.
Uso: python run_quantum.py NOMBRE [NCAS ...]"""
import sys, os, json, time
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(here))
import numpy as np
from pyscf import gto, dft
from psscreen.geometry import read_xyz, to_pyscf_atom
from psscreen.backends import _ecp
from psscreen.quantum_pl import (active_space, eigenstates, threshold_projection, isc_proxy,
                                 taper_report, ftqc_resources, qubit_hamiltonian, EV_NM)
from psscreen.rates import k_isc, MLJParams

ref = json.load(open(os.path.join(here, "references.json")))
name = sys.argv[1]
ncas_list = [int(v) for v in sys.argv[2:]] or [4, 6]
sym, xyz = read_xyz(os.path.join(here, "runs", f"{name}.xyz"))
atom = to_pyscf_atom(sym, xyz)
basis = "def2-svp"
mol = gto.M(atom=atom, basis=basis, ecp=_ecp(atom, basis), verbose=0)
mf = dft.RKS(mol, xc="b3lyp").density_fit(); mf.grids.level = 2; mf.kernel()
out = {"name": name, "cas": {}}
lam_exp = ref[name]["lambda_abs_nm"]; E_exp = EV_NM / lam_exp
for n in ncas_list:
    t0 = time.time()
    asp = active_space(mf, n, n, basis)
    H = qubit_hamiltonian(asp)
    st = eigenstates(asp, H)
    r = {"E_S_eV": st["E_S"][:4].tolist(), "E_T_eV": st["E_T"][:4].tolist(),
         "lambda_S1_nm": EV_NM / st["E_S"][1]}
    # absorción: ventana del paper (700-850 nm), ventana en torno al experimento, ventana en torno a S1(CASCI)
    wins = {"paper_700_850nm": (EV_NM / 850, EV_NM / 700),
            "exp_pm0.15eV": (E_exp - 0.15, E_exp + 0.15),
            "casS1_pm0.25eV": (st["E_S"][1] - 0.25, st["E_S"][1] + 0.25)}
    r["absorption"] = {}
    for k, (lo, hi) in wins.items():
        tt = time.time()
        a = threshold_projection(asp, st, lo, hi, edge_ev=0.1, max_degree=2000)
        a["wall_s"] = time.time() - tt
        r["absorption"][k] = {kk: (float(v) if isinstance(v, (float, np.floating)) else v) for kk, v in a.items()}
        print(name, n, k, f"q={a['weight_quantum']:.4f} exact={a['weight_exact_smooth']:.4f} sharp={a['weight_exact_sharp']:.4f} deg={a["degree"]} edge={a["edge_eV"]:.2f}eV t={a['wall_s']:.0f}s", flush=True)
    r["isc"] = []
    for tr in range(min(3, len(st["T"]))):
        p = isc_proxy(asp, st, times=(20.0, 40.0), triplet=tr)
        dE = st["E_T"][tr] - st["E_S"][1]
        p["dE_ST_eV"] = float(dE)
        p["k_mlj_cas_gap"] = float(k_isc(p["V_proxy_cm"][0], dE, MLJParams()))
        r["isc"].append(p)
        print(name, n, f"T{tr+1}", f"V_proxy={p['V_proxy_cm'][0]:.3f} V_exact={p['V_exact_cm']:.3f} cm-1 dE={dE:+.2f} k_proxy(t=20)={p['k'][0]:.2e}", flush=True)
    r["taper"] = taper_report(asp, H)
    r["wall_s"] = time.time() - t0
    out["cas"][n] = r
    print(name, n, "taper", r["taper"], f"t={r['wall_s']:.0f}s", flush=True)
    json.dump(out, open(os.path.join(here, "runs", f"{name}.quantum.json"), "w"), indent=1, default=float)
# recursos FTQC para un espacio activo del tamaño del paper (11 orbitales)
try:
    asp11 = active_space(mf, 11, 12, basis)
    out["ftqc_cas_12_11"] = ftqc_resources(asp11)
except Exception as e:
    out["ftqc_cas_12_11"] = {"error": repr(e)}
print(name, "FTQC CAS(12,11)", out["ftqc_cas_12_11"], flush=True)
json.dump(out, open(os.path.join(here, "runs", f"{name}.quantum.json"), "w"), indent=1, default=float)
