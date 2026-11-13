## Absorption (S0->S1, vertical, B3LYP/def2-SVP, gas phase, GFN2-xTB geometry)

| Molecule | λ calc (nm) | λ exp (nm) | E calc (eV) | E exp (eV) | Error (eV) | f |
|---|---|---|---|---|---|---|
| BODIPY_core | 507 | 503 | 2.45 | 2.46 | -0.02 | 0.38 |
| BDP | 517 | 503 | 2.40 | 2.46 | -0.07 | 0.51 |
| I2-BDP | 548 | 537 | 2.26 | 2.31 | -0.05 | 0.41 |
| azaBDP_Ph4 | 678 | 650 | 1.83 | 1.91 | -0.08 | 0.65 |

Mean error -0.05 eV, standard deviation 0.02 eV.

Systematic correction with leave-one-out cross-validation (constant shift):

| Molecule | LOO-corrected λ (nm) | λ exp (nm) | LOO error (eV) |
|---|---|---|---|
| BODIPY_core | 494 | 503 | +0.05 |
| BDP | 507 | 503 | -0.02 |
| I2-BDP | 535 | 537 | +0.01 |
| azaBDP_Ph4 | 662 | 650 | -0.03 |

## Triplets, SOC and ISC

| Molecule | T1..T3 (eV) | ΔE(S1-T1) calc / exp | SOC S1-Tn (cm⁻¹) | k_r calc / exp (s⁻¹) | k_ISC calc / exp (s⁻¹) | Φ_T calc / exp | Φ_F calc / exp | proxy Σ|V|² |
|---|---|---|---|---|---|---|---|---|
| BODIPY_core | 1.76, 2.91, 3.10 | 0.71 / — | 0.00, 0.01, 2.52 | 2.2e+08 / — | 2.4e-02 / — | 0.00 / ~0 | 0.60 / ~1 (near unity) | 9.82 |
| BDP | 1.71, 2.91, 3.07 | 0.72 / 0.77 | 0.00, 0.00, 1.67 | 2.9e+08 / 1.2e+08 | 1.1e-02 / not observed (<~1e7) | 0.00 / not observed | 0.66 / 0.45 | 2.84 |
| I2-BDP | 1.68, 2.60, 2.73 | 0.61 / 0.59 | 13.24, 29.39, 8.82 | 2.0e+08 / 1.3e+08 | 6.9e+09 / 3.5e+09 | 0.95 / 0.84 | 0.03 / 0.03 | 1254112.00 |
| azaBDP_Ph4 | 1.15, 2.17, 2.26 | 0.73 / — | 0.11, 0.52, 0.98 | 2.1e+08 / — | 1.5e+05 / — | 0.00 / low (not quantified here) | 0.58 / 0.34 | 5.88 |

### Dominant S1→Tn channels

- BODIPY_core: T1 (ΔE=-0.71 eV, V=0.00 cm⁻¹, k=2.4e-02 s⁻¹); T2 (ΔE=+0.44 eV, V=0.01 cm⁻¹, k=1.0e-08 s⁻¹)
- BDP: T1 (ΔE=-0.72 eV, V=0.00 cm⁻¹, k=1.1e-02 s⁻¹); T2 (ΔE=+0.48 eV, V=0.00 cm⁻¹, k=4.5e-11 s⁻¹)
- I2-BDP: T1 (ΔE=-0.61 eV, V=13.24 cm⁻¹, k=6.9e+09 s⁻¹); T2 (ΔE=+0.32 eV, V=29.39 cm⁻¹, k=1.6e+04 s⁻¹)
- azaBDP_Ph4: T1 (ΔE=-0.73 eV, V=0.11 cm⁻¹, k=1.5e+05 s⁻¹); T2 (ΔE=+0.29 eV, V=0.52 cm⁻¹, k=3.0e+01 s⁻¹)

## Sensitivity of the k_ISC(I2-BDP)/k_ISC(BDP) ratio to the MLJ parameters

Experiment: 3.5e9 vs 'not observed' (Φ_F 0.45 with no detectable ISC, bound k_ISC ≲ 1e7 s⁻¹).

| λ_s (eV) | ħω (eV) | S | k_ISC BDP | k_ISC I2-BDP | Φ_T I2-BDP | ratio |
|---|---|---|---|---|---|---|
| 0.05 | 0.17 | 0.5 | 1.5e-03 | 2.3e+09 | 0.87 | 1.6e+12 |
| 0.05 | 0.17 | 0.8 | 6.2e-03 | 7.1e+09 | 0.95 | 1.1e+12 |
| 0.05 | 0.17 | 1.2 | 1.9e-02 | 1.6e+10 | 0.98 | 8.3e+11 |
| 0.1 | 0.17 | 0.5 | 5.0e-03 | 3.2e+09 | 0.90 | 6.4e+11 |
| 0.1 | 0.17 | 0.8 | 1.6e-02 | 8.0e+09 | 0.96 | 5.1e+11 |
| 0.1 | 0.17 | 1.2 | 3.8e-02 | 1.6e+10 | 0.98 | 4.2e+11 |
| 0.2 | 0.17 | 0.5 | 1.8e-02 | 1.1e+10 | 0.97 | 6.1e+11 |
| 0.2 | 0.17 | 0.8 | 4.0e-02 | 1.9e+10 | 0.98 | 4.8e+11 |
| 0.2 | 0.17 | 1.2 | 7.4e-02 | 2.8e+10 | 0.99 | 3.8e+11 |

## Quantum backend (PennyLane, simulator): CASCI active spaces, no dynamic correlation

| Molecule | CAS | qubits (tapered) | λ(S1) CASCI (nm) | λ exp (nm) | window weight, paper 700-850 nm | window weight around S1 (QSVT / exact) | ISC proxy: best channel V (cm⁻¹), QSVT-Trotter vs exact |
|---|---|---|---|---|---|---|---|
| BODIPY_core | (4,4) | 8 (6) | 266 | 503 | 0.000 | 0.724 / 0.724 | T3: 1.420 vs 1.420 |
| BDP | (4,4) | 8 (6) | 259 | 503 | 0.000 | 0.940 / 0.940 | T3: 0.806 vs 0.806 |
| I2-BDP | (4,4) | 8 (6) | 270 | 537 | 0.000 | 0.847 / 0.847 | T2: 27.914 vs 27.915 |
| I2-BDP | (6,6) | 12 (10) | 272 | 537 | 0.000 | 0.610 / 0.610 | T2: 41.224 vs 41.230 |
| azaBDP_Ph4 | (4,4) | 8 (6) | 336 | 650 | 0.000 | 0.763 / 0.763 | T2: 0.449 vs 0.449 |

ISC proxy ranking (largest S1→Tn V over the first three triplets, CAS(4,4)):

- BODIPY_core: T1 V=0.000 cm⁻¹ (ΔE=-2.70 eV), T2 V=0.008 cm⁻¹ (ΔE=-0.06 eV), T3 V=1.420 cm⁻¹ (ΔE=+1.92 eV)
- BDP: T1 V=0.000 cm⁻¹ (ΔE=-2.79 eV), T2 V=0.004 cm⁻¹ (ΔE=+0.07 eV), T3 V=0.806 cm⁻¹ (ΔE=+2.18 eV)
- I2-BDP: T1 V=11.494 cm⁻¹ (ΔE=-2.36 eV), T2 V=27.915 cm⁻¹ (ΔE=+0.45 eV), T3 V=22.247 cm⁻¹ (ΔE=+2.59 eV)
- azaBDP_Ph4: T1 V=0.067 cm⁻¹ (ΔE=-2.28 eV), T2 V=0.449 cm⁻¹ (ΔE=+0.23 eV), T3 V=0.316 cm⁻¹ (ΔE=+2.44 eV)

FTQC resource estimate for a CAS(12,11) active space (PennyLane double factorization, error 1.6 mHa):

- I2-BDP: 601 logical qubits, 4.39e+07 Toffolis per QPE run, λ = 11.5 Ha
- azaBDP_Ph4: 601 logical qubits, 4.09e+07 Toffolis per QPE run, λ = 10.8 Ha

## Paper-sized active space (12e,11o) on BDP (classical reference)

- Experiment: 503 nm (2.46 eV)
- CASCI (12e,11o): S1 = 4.14 eV (300 nm)
- CASCI + SC-NEVPT2: S1 = 2.52 eV (492 nm)

## Paper-sized active space (12e,11o) on BODIPY_core (classical reference)

- Experiment: 503 nm (2.46 eV)
- CASCI (12e,11o): S1 = 3.97 eV (312 nm)
- CASCI + SC-NEVPT2: S1 = 2.71 eV (458 nm)

## Paper-sized active space (12e,11o) on I2-BDP (classical reference)

- Experiment: 537 nm (2.31 eV)
- CASCI (12e,11o): S1 = 4.32 eV (287 nm)
- CASCI + SC-NEVPT2: S1 = 2.04 eV (609 nm)

## Paper-sized active space (12e,11o) on azaBDP_Ph4 (classical reference)

- Experiment: 650 nm (1.91 eV)
- CASCI (12e,11o): S1 = 3.35 eV (370 nm)
- CASCI + SC-NEVPT2: S1 = 1.99 eV (622 nm)

Wall time per molecule (s): BODIPY_core 241, BDP 2589, I2-BDP 2473, azaBDP_Ph4 3775
