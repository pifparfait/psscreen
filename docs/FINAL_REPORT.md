# psscreen: final report

Status of the framework after validation, a first screening round and a functional scan.
Every number below comes from files in this repository (`validation/`, `screening/`, branch `results-tars`).

## What was built

`psscreen` takes a SMILES string and returns a predicted absorption wavelength, triplet energy, k_ISC, a triplet yield
and a photodynamic-therapy score for BODIPY-type photosensitizers. The route is: RDKit/MMFF → GFN2-xTB geometry →
B3LYP/def2-SVP (PySCF) → own RPA singlets and TDA triplets → ΔSCF correction of S1 and T1 → effective one-electron
spin–orbit coupling → Fermi golden rule with a Marcus–Levich–Jortner Franck–Condon term. Methodology: `docs/METHODOLOGY.md`.

A PennyLane backend implements the approach of arXiv 2512.15889 (QSVT threshold projection, Trotter ISC proxy,
Z2 tapering, double-factorization resource estimates). It runs on simulators only and is not used for rankings.

## What was validated

Nine dyes with experimental data (toluene; Ly et al., PCCP 2021, and others), GFN2-xTB geometries.

| Class | Dyes | Absorption error (calc − exp) | ISC |
|---|---|---|---|
| Compact BODIPY / aza-BODIPY | 6 | −0.02 to −0.08 eV (mean −0.055) | heavy-atom trend reproduced: I2-BDP 6.9e9 vs 3.5e9 s⁻¹, I3-BDP 7.1e9 vs 3.9e9; non-ISC dyes stay at ≈0 |
| Distyryl (BDP-0, I-BDP-0, I2-BDP-0) | 3 | −0.21 to −0.23 eV (λ 77–86 nm too long) | I2-BDP-0: 1.7e7 vs 4.5e8 s⁻¹ (26× low), Φ_T 0.04 vs 0.76 |

**The pipeline is usable for compact dyes and not for extended ones.**

## Why the styryl dyes fail

- Conformers are not the cause (E,E conformers: λ spread 15–18 nm, k_ISC changes <1.5×). An earlier run used Z-styryl
  geometries by mistake; those results were withdrawn (`docs/VALIDATION.md`).
- Absorption: the ΔSCF shift on S1 is roughly constant (≈ −0.5 eV) while the raw TDDFT error shrinks for larger
  chromophores, so the correction overshoots.
- ISC: in I2-BDP-0 the S1→T2 gap is ≈ 0.4–0.5 eV in the corrected scheme and k_ISC through T1 is tiny (SOC ≈ 0.4 cm⁻¹).
  T2 is not shifted like S1 and T1, so the S1–T2 gap is inconsistent. Tests that did not fix it:
  shifting all T_n like S1 (breaks the compact dyes); raw uncorrected TDDFT (reproduces I2-BDP-0 at 1.1e9 s⁻¹ but
  ruins the six compact dyes, so it is cancellation of errors); a hand-tuned T2 offset (works only inside a ~0.1 eV window).

## Functional scan (B3LYP 20 % HF, PBE0 25 %, BHandHLYP 50 %)

| Class | B3LYP | PBE0 | BHandHLYP |
|---|---|---|---|
| 6 compact dyes, absorption error | −0.02 … −0.08 eV | 0.00 … +0.05 eV | +0.34 … +0.38 eV |
| 3 styryl dyes, absorption error | −0.21 … −0.23 eV | −0.14 … −0.15 eV | +0.21 … +0.24 eV |

- PBE0 improves absorption for both classes. Linear interpolation puts the zero error near 24 % HF (compact) and
  near 35 % HF (styryl), which hints at a size-dependent optimal exchange fraction. This rests on three functionals
  and nine dyes and was not tested at intermediate fractions.
- PBE0 keeps the iodinated compact dyes reasonable (I2-BDP 1.6e9 vs 3.5e9; I3-BDP 1.6e9 vs 3.9e9). BHandHLYP does not (~2e6).
- The styryl ISC failure does not depend on the functional: I2-BDP-0 gives 1.7e7 (B3LYP), 7.6e6 (PBE0), 1.5e7 (BHandHLYP) vs 4.5e8.
- Only global hybrids are supported by the in-house TDDFT code; range-separated functionals (CAM-B3LYP, ωB97X-D) were not tested.

## Screening round 1 (predictions, no experimental data)

| Candidate | λ (nm) | Φ_T | Comment |
|---|---|---|---|
| azaBDP_Br2 | 687 | 0.97 | **Fragile**: Φ_T comes from an S1→T2 channel with T2 only 0.11 eV above S1 and uncorrected TDDFT energies; through T1 alone it would rank last of the three aza dyes |
| azaBDP_IPh2 | 695 | 0.54 | T1 channel, SOC 5 cm⁻¹ |
| azaBDP_BrPh2 | 701 | 0.04 | |
| Cl2-/Br2-/I1-/I3-BDP | ≈550 | ≈0 PDT score | outside the 650–850 nm window; k_ISC trend Cl2 < Br2 < I1 < I3 as expected |

Aza candidates may be somewhat too red (validation error −0.08 eV on the aza reference). Treat the ranking as a
shortlist for a more reliable method, not a result.

## Quantum backend, honest assessment

CASCI(12,11) puts S1 1.4–1.9 eV too high; NEVPT2 is within 0.06–0.27 eV with sign changes; the Trotter ISC proxy has no
S–T gap. A fault-tolerant estimate for CAS(12,11) is ~601 logical qubits and ~4e7 Toffolis. No new method: the
contribution is a benchmark (negative for the paper's observables as screening descriptors), a validated classical
pipeline and an open framework.

## Limits

- Valid only for compact BODIPY / aza-BODIPY dyes; extended π-systems fail.
- k_ISC depends on S1–Tn gaps at the 0.1 eV level, comparable to the method's accuracy; heavy-atom ISC channels via T2 are not predictive.
- Single conformer per dye except the distyryl set; gas-phase TDDFT with fixed k_IC (1.5e8 s⁻¹).
- Not tested: macOS Apple Silicon install path.

## Reproduce

`scripts/run_validation_server.sh`, `scripts/run_screen_parallel.sh`, `scripts/run_xc_scan_server.sh` (Linux server);
`README.md` for a local run. Results of server runs are on branch `results-tars`.

## Suggested next steps (not done)

1. Treat T1-based quantities as the validated output; report S1→T2 channels only with an explicit gap uncertainty.
2. Test 30 and 35 % HF and a range-separated functional (requires adding attenuated exchange to `tdwin.py`).
3. A size-dependent S1 correction, fitted on more extended dyes than the three used here.
4. Experimental check of the top aza candidates before any further computation.
