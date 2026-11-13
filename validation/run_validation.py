"""Ejecuta el pipeline sobre el conjunto de validación. Uso: python run_validation.py [nombres...]"""
import sys, os, json
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(here)); sys.path.insert(0, here)
from psscreen.pipeline import screen_molecule
from psscreen.backends import PySCFTDA
from molecules import SMILES

names = sys.argv[1:] or list(SMILES)
for n in names:
    r = screen_molecule(n, SMILES[n], backend=PySCFTDA(xc="b3lyp", basis="def2-svp", nstates=5),
                        workdir=os.path.join(here, "runs"), n_refr=1.496)  # tolueno
    p = r["photophysics"]
    print(f"{n}: lam={r['lambda_max_nm']:.0f} nm f={r['f_osc'][0]:.2f} S1={r['E_S1_eV']:.2f} T={['%.2f'%e for e in r['E_triplets_eV']]} "
          f"SOC={['%.2f'%v for v in r['SOC_S1_Tn_cm']]} kISC={p["k_isc"]:.2e} raw_lam={r["raw_tddft"]["lambda_nm"]:.0f} raw_phiISC={r["raw_tddft"]["photophysics"]["phi_ISC"]:.2f} phiISC={p['phi_ISC']:.2f} t={r['wall_s']['total']:.0f}s", flush=True)
