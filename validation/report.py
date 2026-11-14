"""Build the predicted-vs-experimental tables and the sensitivity analysis from runs/*.json."""
import json, os, sys, itertools
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(here))
from psscreen.rates import photophysics, MLJParams

ref = json.load(open(os.path.join(here, "references.json")))
def _path(n):
    """Result file of a dye; flexible dyes use their lowest-energy conformer (set3/<name>_c0.json)."""
    for cand in (f"{n}.json", os.path.join("set3", f"{n}_c0.json")):
        if os.path.exists(os.path.join(here, "runs", cand)):
            return os.path.join(here, "runs", cand)
    return None


names = [n for n in ref if _path(n)]
R = {n: json.load(open(_path(n))) for n in names}
EV = 1239.84198
rows = []
for n in names:
    r, e = R[n], ref[n]; p = r["photophysics"]
    lam_exp = e["lambda_abs_nm"]; E_exp = EV / lam_exp
    rows.append(dict(name=n, lam=r["lambda_max_nm"], lam_exp=lam_exp, E=r["E_S1_eV"], E_exp=E_exp,
                     dE=r["E_S1_eV"] - E_exp, f=r["f_osc"][0], ET=r["E_triplets_eV"], dEST=r["dE_ST_eV"],
                     dEST_exp=e.get("dE_ST_eV"), soc=r["SOC_S1_Tn_cm"], kisc=p["k_isc"], kisc_exp=e.get("k_isc_s"),
                     phiT=p["phi_ISC"], phiT_exp=e.get("phi_T"), phiF=p["phi_F"], phiF_exp=e.get("phi_F"),
                     kr=p["k_r"], kr_exp=e.get("k_r_s"), proxy=r["proxy_xanadu_sum_V2"], ch=p["channels"], t=r["wall_s"]["total"]))

print("## Absorption (S0->S1, vertical, B3LYP/def2-SVP, gas phase, GFN2-xTB geometry)\n")
print("| Molecule | λ calc (nm) | λ exp (nm) | E calc (eV) | E exp (eV) | Error (eV) | f |")
print("|---|---|---|---|---|---|---|")
for x in rows:
    print(f"| {x['name']} | {x['lam']:.0f} | {x['lam_exp']} | {x['E']:.2f} | {x['E_exp']:.2f} | {x['dE']:+.2f} | {x['f']:.2f} |")
errs = np.array([x["dE"] for x in rows])
print(f"\nMean error {errs.mean():+.2f} eV, standard deviation {errs.std():.2f} eV.")
if len(rows) >= 3:
    print("\nSystematic correction with leave-one-out cross-validation (constant shift):\n")
    print("| Molecule | LOO-corrected λ (nm) | λ exp (nm) | LOO error (eV) |\n|---|---|---|---|")
    for i, x in enumerate(rows):
        sh = np.mean([y["dE"] for j, y in enumerate(rows) if j != i])
        Ec = x["E"] - sh
        print(f"| {x['name']} | {EV / Ec:.0f} | {x['lam_exp']} | {Ec - x['E_exp']:+.2f} |")

print("\n## Triplets, SOC and ISC\n")
print("| Molecule | T1..T3 (eV) | ΔE(S1-T1) calc / exp | SOC S1-Tn (cm⁻¹) | k_r calc / exp (s⁻¹) | k_ISC calc / exp (s⁻¹) | Φ_T calc / exp | Φ_F calc / exp | proxy Σ|V|² |")
print("|---|---|---|---|---|---|---|---|---|")
fmt = lambda v: "—" if v is None else (f"{v:.1e}" if isinstance(v, float) and (v > 100 or v < 1e-2) else str(v))
for x in rows:
    print(f"| {x['name']} | {', '.join(f'{e:.2f}' for e in x['ET'][:3])} | {x['dEST']:.2f} / {fmt(x['dEST_exp'])} | "
          f"{', '.join(f'{v:.2f}' for v in x['soc'][:3])} | {x['kr']:.1e} / {fmt(x['kr_exp'])} | {x['kisc']:.1e} / {fmt(x['kisc_exp'])} | "
          f"{x['phiT']:.2f} / {x['phiT_exp']} | {x['phiF']:.2f} / {x['phiF_exp']} | {x['proxy']:.2f} |")

print("\n### Dominant S1→Tn channels\n")
for x in rows:
    top = sorted(x["ch"], key=lambda c: -c["k"])[:2]
    print(f"- {x['name']}: " + "; ".join(f"T{c['T']} (ΔE={c['dE_eV']:+.2f} eV, V={c['V_cm']:.2f} cm⁻¹, k={c['k']:.1e} s⁻¹)" for c in top))

if "BDP" in R and "I2-BDP" in R:
    print("\n## Sensitivity of the k_ISC(I2-BDP)/k_ISC(BDP) ratio to the MLJ parameters\n")
    print("Experiment: 3.5e9 vs 'not observed' (Φ_F 0.45 with no detectable ISC, bound k_ISC ≲ 1e7 s⁻¹).\n")
    print("| λ_s (eV) | ħω (eV) | S | k_ISC BDP | k_ISC I2-BDP | Φ_T I2-BDP | ratio |\n|---|---|---|---|---|---|---|")
    for ls, hw, S in itertools.product([0.05, 0.1, 0.2], [0.17], [0.5, 0.8, 1.2]):
        p = MLJParams(lambda_s=ls, hw=hw, S=S)
        a = photophysics(R["BDP"]["E_S1_eV"], R["BDP"]["f_osc"][0], R["BDP"]["E_triplets_eV"], R["BDP"]["SOC_S1_Tn_cm"], p=p, n_refr=1.496)
        b = photophysics(R["I2-BDP"]["E_S1_eV"], R["I2-BDP"]["f_osc"][0], R["I2-BDP"]["E_triplets_eV"], R["I2-BDP"]["SOC_S1_Tn_cm"], p=p, n_refr=1.496)
        print(f"| {ls} | {hw} | {S} | {a['k_isc']:.1e} | {b['k_isc']:.1e} | {b['phi_ISC']:.2f} | {b['k_isc'] / max(a['k_isc'], 1e-30):.1e} |")

# ---------------------------------------------------------------------------
# Quantum backend (PennyLane simulator) vs classical vs experiment
qnames = [n for n in ["BODIPY_core", "BDP", "I2-BDP", "azaBDP_Ph4"]
          if os.path.exists(os.path.join(here, "runs", f"{n}.quantum.json"))]
if qnames:
    Q = {n: json.load(open(os.path.join(here, "runs", f"{n}.quantum.json"))) for n in qnames}
    print("\n## Quantum backend (PennyLane, simulator): CASCI active spaces, no dynamic correlation\n")
    print("| Molecule | CAS | qubits (tapered) | λ(S1) CASCI (nm) | λ exp (nm) | window weight, paper 700-850 nm | window weight around S1 (QSVT / exact) | ISC proxy: best channel V (cm⁻¹), QSVT-Trotter vs exact |")
    print("|---|---|---|---|---|---|---|---|")
    for n in qnames:
        for c, r in Q[n]["cas"].items():
            a = r["absorption"]
            w700 = a["paper_700_850nm"]["weight_quantum"]
            ws = a["casS1_pm0.25eV"]
            best = max(r["isc"], key=lambda p: p["V_exact_cm"])
            print(f"| {n} | ({c},{c}) | {r['taper']['qubits']} ({r['taper']['qubits_tapered']}) | {r['lambda_S1_nm']:.0f} | "
                  f"{ref[n]['lambda_abs_nm']} | {w700:.3f} | {ws['weight_quantum']:.3f} / {ws['weight_exact_smooth']:.3f} | "
                  f"T{best['triplet']}: {best['V_proxy_cm'][0]:.3f} vs {best['V_exact_cm']:.3f} |")
    print("\nISC proxy ranking (largest S1→Tn V over the first three triplets, CAS(4,4)):\n")
    for n in qnames:
        r = Q[n]["cas"].get("4")
        if r:
            print(f"- {n}: " + ", ".join(f"T{p['triplet']} V={p['V_exact_cm']:.3f} cm⁻¹ (ΔE={p['dE_ST_eV']:+.2f} eV)" for p in r["isc"]))
    fq = [(n, Q[n]["ftqc_cas_12_11"]) for n in qnames if "logical_qubits" in Q[n].get("ftqc_cas_12_11", {})]
    if fq:
        print("\nFTQC resource estimate for a CAS(12,11) active space (PennyLane double factorization, error 1.6 mHa):\n")
        for n, f in fq:
            print(f"- {n}: {f['logical_qubits']} logical qubits, {f['toffoli_gates']:.2e} Toffolis per QPE run, λ = {f['lambda']:.1f} Ha")


# Paper-sized active space reference (classical CASCI / NEVPT2)
import glob
for f in sorted(glob.glob(os.path.join(here, "runs", "*.casci_ref.json"))):
    c = json.load(open(f))
    print(f"\n## Paper-sized active space {c['CAS']} on {c['name']} (classical reference)\n")
    print(f"- Experiment: {c['lambda_exp_nm']} nm ({1239.84198 / c['lambda_exp_nm']:.2f} eV)")
    print(f"- CASCI {c['CAS']}: S1 = {c['E_CASCI_eV'][1]:.2f} eV ({c['lambda_CASCI_nm'][0]:.0f} nm)")
    if "E_NEVPT2_eV" in c:
        print(f"- CASCI + SC-NEVPT2: S1 = {c['E_NEVPT2_eV'][1]:.2f} eV ({c['lambda_NEVPT2_nm'][0]:.0f} nm)")

print("\nWall time per molecule (s): " + ", ".join(f"{x['name']} {x['t']:.0f}" for x in rows))
