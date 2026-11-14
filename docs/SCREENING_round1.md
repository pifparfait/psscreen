# Screening round 1: heavy-atom and red-shifted BODIPY / aza-BODIPY photosensitizers

Predictions from the validated pipeline (GFN2-xTB geometry, B3LYP/def2-SVP, own RPA singlets, ΔSCF
S1/T1 correction, one-electron SOC, Marcus-Levich-Jortner k_ISC; see `docs/VALIDATION.md`).
Window for the score: 650-850 nm. Candidates are in `screening/candidates_round1.csv`; per-molecule
JSON and geometries are in `screening/runs/candidates_round1/`.

> **Correction (2026-11-14).** The first version of this report listed an "I2-distyryl-BDP" as the top
> candidate (673 nm, k_ISC 2.1e10 s⁻¹, Φ_T 0.98). Its SMILES was wrong: it put one styryl group at
> position 1 and the other at position 5 instead of 3,5. That result was discarded. The correct molecule
> (3,5-distyryl, 2,6-diiodo, meso-phenyl) is "I2-BDP-0" in Ly et al., PCCP 2021, which reports 644 nm,
> k_ISC = 4.5e8 s⁻¹ and Φ_T = 0.76 in toluene. It is now part of the second validation set instead
> of the candidate list, and the results below cover the remaining seven candidates.

## Results

| Candidate | λ_max (nm) | f | S1 (eV) | T1 (eV) | ΔE(S1-T1) (eV) | k_ISC (s⁻¹) | Φ_T | PDT score |
|---|---|---|---|---|---|---|---|---|
| azaBDP_Br2 | 687 | 0.66 | 1.80 | 1.12 | 0.68 | 1.3e10 | 0.97 | 0.97 |
| azaBDP_IPh2 | 695 | 0.66 | 1.78 | 1.08 | 0.71 | 4.2e8 | 0.54 | 0.54 |
| azaBDP_BrPh2 | 701 | 0.67 | 1.77 | 1.07 | 0.70 | 1.5e7 | 0.04 | 0.04 |
| I3-BDP | 551 | 0.40 | 2.25 | 1.64 | 0.61 | 7.1e9 | 0.95 | 0.02 |
| Br2-BDP | 548 | 0.40 | 2.26 | 1.65 | 0.61 | 3.3e8 | 0.48 | 0.01 |
| I1-BDP | 531 | 0.38 | 2.34 | 1.67 | 0.67 | 2.3e9 | 0.87 | 0.00 |
| Cl2-BDP | 550 | 0.37 | 2.26 | 1.63 | 0.62 | 2.6e6 | 0.01 | 0.00 |

For reference, the validated I2-BDP gave k_ISC = 6.9e9 s⁻¹ (experiment 3.5e9), Φ_T = 0.95.

## What the numbers say

- **Heavy-atom trend inside the BDP scaffold is physically consistent**: k_ISC rises Cl2 (2.6e6) < Br2 (3.3e8)
  < I1 (2.3e9) < I3 (7.1e9) s⁻¹. These BDP variants all absorb near 550 nm, outside the 650-850 nm
  window, so they score about 0 and serve as trend checks, not candidates.
- **Only the 2,6-dibromo aza-BODIPY combines absorption in the window (687 nm) with a high predicted
  triplet yield (0.97).** The distyryl-diiodo dye that first led the ranking was removed (see correction above).
- **Where the heavy atom sits matters**: bromine on the pendant phenyls (azaBDP_BrPh2) gives almost no ISC
  (Φ_T 0.04), while iodine on the same positions gives a moderate yield (0.54) and bromine on the
  pyrrole core gives 0.97.

## Caveats (read before acting on the ranking)

1. **These are predictions, not results.** Only four molecules were validated, and only one (I2-BDP) has
   a measured triplet yield. No candidate here has experimental data in this work.
2. **Systematic λ error.** Over nine validation dyes the calculated absorption is on average 0.07 eV too
   long (+4 to +50 nm). The aza candidates at 687-701 nm may realistically lie around 660-675 nm; they
   should still fall inside the window.
3. **Styryl-extended dyes are not reliable yet.** The second validation set shows that for the
   distyryl-diiodo dye the method gets Φ_T = 0.18 against 0.76 measured, and that the absorption of
   styryl dyes depends strongly on the conformer (see `docs/VALIDATION.md`). No styryl candidate is in
   this round for that reason.
4. **T1 margin for singlet oxygen.** The aza-BODIPY triplets sit at 1.07-1.12 eV, only 0.09-0.14 eV above
   the 0.98 eV of O2(1Δg). The ΔSCF error is probably of that size, so type II sensitization is not
   guaranteed for these candidates. The score only checks T1 > 1.03 eV, so it does not reflect this.
5. **k_ISC depends on empirical parameters.** The reorganization and vibrational parameters were set a
   priori; for I2-BDP they move k_ISC between about 2e9 and 3e10 s⁻¹. The ordering Cl < Br < I is
   robust to that, but the exact Φ_T of Br2-BDP (0.48) and azaBDP_IPh2 (0.54) is not: both sit where
   k_ISC is comparable to k_r + k_ic.
6. **k_ic is fixed at 1.5e8 s⁻¹** (taken from the BDP data), so Φ_F and Φ_T of the new dyes inherit that choice.
7. **Gas-phase, single conformer, GFN2-xTB geometry**, no solvent or aggregation effects.

## Suggested next steps

- The second validation set is done (see `docs/VALIDATION.md`). Next: test whether applying the ΔSCF
  correction (or a better treatment) to the higher triplets fixes the distyryl-diiodo k_ISC, add conformer
  sampling for styryl dyes, and test a range-separated functional and a solvent model.
- Look for measured data on 2,6-dihalogenated aza-BODIPYs to validate the top remaining candidate.
- Re-run the top candidates with an optimized (DFT) geometry.
