# Methodology

What `psscreen` computes, step by step, and what is known to be approximate. All steps are in `psscreen/`.

## Classical pipeline (produces the rankings)

1. **Geometry** (`geometry.py`). SMILES → RDKit embedding and MMFF → GFN2-xTB optimization (tblite, ASE).
   One conformer, gas phase.
2. **Ground state and excited states** (`backends.py`, `tdwin.py`). B3LYP/def2-SVP with PySCF and density
   fitting; an effective core potential only for iodine.
   - Singlets: an own RPA (full TDDFT) Davidson solver in the MO basis, using the Stratmann-Scuseria-Frisch
     scheme, restricted to an orbital window (up to a 60 eV energy window).
   - Triplets: explicit TDA in an orbital window (`window_td`), validated against PySCF.
3. **ΔSCF correction** (`deltascf.py`). Raw TDDFT puts S1 0.35-0.76 eV too high for the cyanine-like
   BODIPY chromophore and produces spurious triplet yields. T1 comes from an unrestricted M_S=1 calculation.
   S1 comes from a mixed determinant with maximum-overlap-method (MOM) convergence and the Ziegler sum
   rule, E(S) = 2 E(mixed) − E(T).
4. **Spin-orbit coupling** (`soc.py`). Effective one-electron spin-orbit operator with a per-element scale
   factor calibrated so that the isolated atom reproduces the experimental ζ(np). Matrix elements between
   S1 and each triplet (TDA states) are summed over M_S. Checked against the El-Sayed rules in formaldehyde.
5. **Rates** (`rates.py`).
   - k_ISC: Fermi golden rule with a Marcus-Levich-Jortner Franck-Condon weighted density of states
     (λ_s = 0.10 eV, ħω = 0.17 eV, S = 0.8), summed over the S–T gap and all triplets.
   - k_r: Einstein / Strickler-Berg from the computed oscillator strength and the solvent refractive index.
   - k_ic: fixed at 1.5e8 s⁻¹ (taken from the measured non-radiative rate of BDP).
   - Φ_T = k_ISC / (k_ISC + k_r + k_ic).
6. **Score** (`pipeline.py`, `pdt_score`). Absorption in the 650-850 nm window × oscillator-strength
   factor × Φ_T × (T1 above the O2 singlet energy, 0.98 eV + 0.05 eV).

## Quantum backend (`quantum_pl.py`, PennyLane simulator)

Follows arXiv:2512.15889: CAS(n,n) → 2n qubits (Jordan-Wigner); threshold projection with QSVT
(`qml.QSVT`, `PCPhase`) and a block encoding; ISC proxy with `qml.TrotterProduct`; Z2 tapering; double
factorization resource estimates for CAS(12,11). In simulation these reproduce exact values. They are not
used for the rankings: CASCI(12e,11o) puts S1 1.4-1.9 eV too high, and CASCI + SC-NEVPT2 is within
0.06-0.27 eV with a sign that changes between molecules (`docs/VALIDATION.md`).

## What is validated, and what is not

- Validated against measurements: absorption of four dyes (errors of +4 to +28 nm), T1 and ΔE(S1-T1) for the BDP series,
  one measured k_ISC (I2-BDP: 6.9e9 predicted vs 3.5e9 measured).
- Not independent: Φ_F (k_ic is taken from BDP).
- Assumed, not fitted: the MLJ parameters. For I2-BDP they move k_ISC between about 2e9 and 3e10 s⁻¹.
- Not included: solvent, conformers, vibronic structure of the absorption band, aggregation.
- A second validation set from the same source (toluene) is in `docs/VALIDATION.md`: the absorption offset
  is systematic (mean -0.055 eV over six dyes). The three styryl-extended dyes are being recalculated after a
  geometry error (missing E/Z labels in the SMILES) invalidated a first attempt.
- Still to test: range-separated functionals, a solvent model, conformer sampling, and ΔSCF-quality energies
  for higher triplets.
