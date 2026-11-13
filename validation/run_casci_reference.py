"""Classical CASCI(12e,11o) reference for a paper-sized active space (no quantum hardware needed).

Tests the main open question about arXiv:2512.15889: how accurate is the S1 energy of an
active space of that size, without dynamic correlation (and, optionally, with SC-NEVPT2)?
Usage: python run_casci_reference.py NAME [--nevpt2]
"""
import sys, os, json, time
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(here))
from pyscf import gto, scf, mcscf, fci
from psscreen.geometry import read_xyz, to_pyscf_atom
from psscreen.backends import _ecp

name = sys.argv[1]; do_nevpt2 = "--nevpt2" in sys.argv
ref = json.load(open(os.path.join(here, "references.json")))
sym, xyz = read_xyz(os.path.join(here, "runs", f"{name}.xyz"))
atom = to_pyscf_atom(sym, xyz)
mol = gto.M(atom=atom, basis="def2-svp", ecp=_ecp(atom, "def2-svp"), verbose=0,
            max_memory=float(os.environ.get("PSS_MAX_MEMORY_MB", 4000)))
t0 = time.time()
mf = scf.RHF(mol).density_fit().run()
mc = mcscf.CASCI(mf, 11, 12)
mc.fcisolver = fci.direct_spin0.FCI(mol); mc.fcisolver.nroots = 4
mc.kernel()
EV = 27.211386245988
E = [(e - mc.e_tot[0]) * EV for e in mc.e_tot]
out = {"name": name, "CAS": "(12e,11o)", "E_CASCI_eV": E, "lambda_CASCI_nm": [1239.84198 / e for e in E[1:]],
       "lambda_exp_nm": ref[name]["lambda_abs_nm"], "t_s": time.time() - t0}
print(name, "CASCI(12,11) excitation energies (eV):", [round(e, 2) for e in E], flush=True)
if do_nevpt2:
    from pyscf import mrpt
    try:
        e2 = [mc.e_tot[i] + mrpt.NEVPT(mc, root=i).kernel() for i in range(3)]
        E2 = [(e - e2[0]) * EV for e in e2]
        out["E_NEVPT2_eV"] = E2; out["lambda_NEVPT2_nm"] = [1239.84198 / e for e in E2[1:]]
        print(name, "NEVPT2 excitation energies (eV):", [round(e, 2) for e in E2], flush=True)
    except RuntimeError as exc:  # keep the CASCI result if NEVPT2 runs out of memory
        out["nevpt2_error"] = str(exc)
        print(name, "NEVPT2 failed:", exc, flush=True)
json.dump(out, open(os.path.join(here, "runs", f"{name}.casci_ref.json"), "w"), indent=1)
print(name, "done", round(time.time() - t0), "s", flush=True)
