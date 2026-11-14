"""Test: does shifting the higher triplets (T2...) like S1 instead of like T1 fix the I2-BDP-0 k_ISC?

Compares the default scheme (all triplets shifted rigidly by the Delta-SCF T1 shift) with a scheme where
T2 and above get the Delta-SCF S1 shift. Uses the result JSONs in runs/, no new quantum chemistry.
Usage: python analysis_higher_triplets.py
"""
import json, os, sys
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(here))
from psscreen.rates import photophysics

NAMES = ["BODIPY_core", "BDP", "I-BDP", "I2-BDP", "I3-BDP", "azaBDP_Ph4", "BDP-0", "I-BDP-0", "I2-BDP-0"]
print("molecule | rigid shift: k_ISC, Phi_T | T2+ shifted like S1: k_ISC, Phi_T")
for n in NAMES:
    d = json.load(open(os.path.join(here, "runs", f"{n}.json"))); c = d["corrections"]
    E = list(d["E_triplets_eV"])
    kw = dict(n_refr=1.496, shift_S=c["S"], shift_T=c["T"])
    a = photophysics(d["E_singlets_eV"][0], d["f_osc"][0], E, d["SOC_S1_Tn_cm"], **kw)
    E2 = [E[0]] + [e + (c["S"] - c["T"]) for e in E[1:]]
    b = photophysics(d["E_singlets_eV"][0], d["f_osc"][0], E2, d["SOC_S1_Tn_cm"], **kw)
    print(f"{n:12s} | {a['k_isc']:.1e} {a['phi_ISC']:.2f} | {b['k_isc']:.1e} {b['phi_ISC']:.2f}")
