# Gap-resolved ISC observable: simulation results (item 2 of the paper plan)

Code: `psscreen/isc_filtered.py`; runner: `paper/quantum/run_filtered_isc.py`; raw results:
`paper/quantum/results/*.filtered_isc.json`. B3LYP/def2-SVP orbitals, CASCI active spaces, default MLJ
parameters (λ_s = 0.10 eV). The QSVT circuit runs on the PennyLane `default.qubit` simulator with an explicit
block encoding in a randomly rotated basis of each M_S sector. "poly" is the same polynomial evaluated on the
exact spectrum, and "exact" is the golden rule with the exact MLJ line shape.

| Dye | CAS | wires (M_S=0 / ±1) | degree | k QSVT (s⁻¹) | k poly | k exact | gap-free S1-T1 proxy V (cm⁻¹) | CASCI ΔE(S1-T1) (eV) |
|---|---|---|---|---|---|---|---|---|
| BODIPY core | (4,4) | 8 / 6 | 3324 | 1.972e4 | 1.972e4 | 1.972e4 | 0.000 | -2.70 |
| BODIPY core | (6,6) | 11 / 10 | 5720 | 2.479e8 | 2.479e8 | 2.479e8 | 0.001 | -2.34 |
| BDP | (4,4) | 8 / 6 | 4106 | 3.519e2 | 3.519e2 | 3.519e2 | 0.000 | -2.79 |
| I2-BDP | (4,4) | 8 / 6 | 4492 | 4.568e-1 | 4.568e-1 | 4.680e-1 | 11.49 | -2.36 |
| I2-BDP | (6,6) | 11 / 10 | 8000 | 4.229e6 | 4.229e6 | 4.230e6 | 24.84 | -2.49 |

## Findings

1. **The circuit computes the observable it is designed for.** The QSVT estimate agrees with the polynomial on
   the exact spectrum to all printed digits, and with the exact golden rule to within 0.03% (2% for I2-BDP
   CAS(4,4)). That rate is 10⁻¹⁰ of the filter maximum, near the dynamic-range floor set by the polynomial
   error (about err² of the maximum rate).
2. **With CASCI energies, the gap-resolved rate ranks the dyes in the wrong order.** In these active spaces S1
   lies 2.3-2.8 eV above T1 instead of about 0.6 eV, so the S1→T1 channel is closed by the Franck-Condon factor.
   The computed rate then comes from upper triplets that happen to lie within 0.06-0.25 eV of the (wrong) S1.
   At CAS(6,6) this puts the iodine-free BODIPY core (2.5×10⁸ s⁻¹) above I2-BDP (4.2×10⁶ s⁻¹), the opposite of
   experiment.
3. **The gap-free proxy ranks these dyes correctly for the wrong reason.** It ignores the energies, so it only
   sees the coupling, which iodine increases (11-25 cm⁻¹ against ≤0.001 cm⁻¹). This is robust to energy errors
   but blind to the mechanism. It fails when an upper triplet near resonance controls ISC (azaBDP_Br2 in the
   screening, see `paper/analysis/accuracy_budget_notes.md`).
4. **Consequence for the paper.** A mechanism-correct ISC observable needs gaps to about 0.07 eV (the accuracy
   budget), and the bare active-space Hamiltonian misses them by about 2 eV. The useful quantum target is
   therefore energies with dynamic correlation at the budget accuracy, not a larger bare active space or a
   better proxy. The filter is the right observable once the Hamiltonian is good enough. It has the same cost
   structure as the threshold projection: d ≈ 4.76 λ'/Δ + 321 with Δ ≈ 0.07 eV, about 3×10⁴ for λ' ≈ 16 Ha.
