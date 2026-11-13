"""Backend PennyLane en formaldehído CAS(4,4): energías = CASCI de PySCF, QSVT = peso exacto,
proxy ISC = |<T|H_SOC|S>|^2 t^2."""
import numpy as np, sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pyscf import gto, dft, mcscf
from psscreen.quantum_pl import active_space, eigenstates, threshold_projection, isc_proxy, taper_report

def test_formaldehyde():
    mol = gto.M(atom="C 0 0 0; O 0 0 1.205; H 0 0.94 -0.587; H 0 -0.94 -0.587", basis="def2-svp", verbose=0)
    mf = dft.RKS(mol, xc="b3lyp").run()
    asp = active_space(mf, 4, 4, "def2-svp")
    st = eigenstates(asp)
    mc = mcscf.CASCI(mf, 4, 4); mc.fcisolver.nroots = 3; mc.kernel()
    print("E_S (eV):", st["E_S"][:3], "E_T:", st["E_T"][:2])
    assert abs(st["E0"] - mc.e_tot[0]) < 1e-7
    lo, hi = st["E_S"][2] - 0.3, st["E_S"][2] + 0.3   # S2 (pi->pi*, permitido)
    r = threshold_projection(asp, st, lo, hi, edge_ev=0.05, max_degree=1500)
    print(r)
    assert abs(r["weight_quantum"] - r["weight_exact_smooth"]) < 1e-4   # error de interpolación del polinomio
    p0 = isc_proxy(asp, st, times=(10.0,), triplet=0)   # S1(n pi*)-T1(n pi*): El-Sayed prohibido
    assert p0["V_exact_cm"] < 1e-6
    p = isc_proxy(asp, st, times=(10.0, 20.0), triplet=1)   # S1(n pi*)-T2(pi pi*): permitido
    print(p)
    print("proxy check")
    assert abs(p["V_proxy_cm"][0] - p["V_exact_cm"]) / max(p["V_exact_cm"], 1e-6) < 1e-2
    print(taper_report(asp))

if __name__ == "__main__":
    test_formaldehyde(); print("OK")
