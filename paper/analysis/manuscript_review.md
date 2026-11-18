# Referee report: "Accuracy Budget, Gap-Resolved Intersystem-Crossing Observable, and Decision-Level Quantum Resource Estimates" (JCTC draft)

Manuscript reviewed: `paper/main.tex` at commit c36a1f7 (branch `claude/project-thread-l7txl1`). Line numbers below
refer to that file. Every number in the abstract, the text and the tables was checked against the generated sources
listed at the end. A number that is not mentioned matched its source.

## Summary and recommendation

The paper makes four contributions. It derives an accuracy budget for an ISC-based screening decision. It
benchmarks a TDDFT/ΔSCF pipeline against measured BODIPY photophysics. It implements and verifies a QSVT
gap-resolved ISC observable. And it costs the screening decision, not the single shot, on a fault-tolerant
machine. The accuracy-budget idea and the decision-level cost model are useful and clearly presented. The
quantum simulation and the cost tables are reproducible from the repository and match their sources almost
everywhere.

My recommendation is **major revision**. Three problems weaken the main thesis as it now stands:

1. The multireference results, which carry the claim that dynamic correlation and not the active space is the
   bottleneck, were run with uncontrolled active spaces.
2. The I2-BDP-0 discussion and the abstract say the methods disagree on the position of T2. The data in the
   paper say the opposite.
3. Several headline ranges in the abstract and the conclusions do not match Table `tab:decision`.

None of these needs new theory. Items 1 and 2 need a few more calculations or a rewrite. Item 3 needs
rewording only.

---

## Major points (ordered by importance)

### M1. The active spaces are not controlled, so the CAS and NEVPT2 conclusions are not yet supported

**Where:** Methods, "Active-space references and resources" (l. 266–276); sec:cas (l. 428–469); Tables
`tab:cas`, `tab:sanevpt2`, `tab:sweep`; abstract l. 54–55; conclusions l. 715–716.

**Problem.** Both `references/run_sa_nevpt2.py` and `references/run_cas_sweep.py` take the active orbitals
straight from the canonical RHF orbitals around the Fermi level. Neither uses AVAS, π-orbital selection or a
check of the orbital character. Three consequences follow.

- **Different orbitals across dyes.** The "(12e,11o)" space then contains different kinds of orbitals in
  different dyes. In the iodinated dyes the iodine lone pairs lie in the HOMO region and can displace BODIPY π
  orbitals. The same label does not describe the same model chemistry across Table `tab:sanevpt2`, and it may
  not be the π space of ref. Zhou2026PRA. The NEVPT2 corrections on SA-CASSCF orbitals point to this. They range
  from +0.11 eV (aza-BDP-Ph4) to −0.80 eV (I-BDP-0), with CASSCF S1 values that do not follow the expected
  trend (BODIPY core 3.43 eV, BDP 3.11 eV).
- **Results depend strongly on the orbitals.** NEVPT2 is 0.06–0.27 eV from experiment on RHF orbitals
  (`tab:cas`) but 0.51–0.79 eV from experiment on SA-CASSCF orbitals (`tab:sanevpt2`), for the same nominal
  active space and the same dyes. For BDP the values are 2.52 against 3.05 eV. SC-NEVPT2 should not move by 0.5
  eV on an orbital change unless the active space changed in character. The paper reports both numbers without
  comment. Instead it attributes the SA error to "a limitation of the second-order treatment" (l. 464–466).
- **The (14,14) step proves little.** The (14,14) step changes every CASSCF/CASCI energy by less than 0.003 eV
  for both dyes (`sweep_partial_2026-11-17.md`). This includes T2–S1 for I2-BDP, which had just jumped by 0.69
  eV between (10,10) and (12,12). The near-identity is the signature of two added orbitals with occupations
  ≈ 2 and ≈ 0. A typical case is a σ or iodine lone pair plus a high σ* orbital that the warm-started SA-CASSCF
  does not rotate. It is not evidence of convergence. "The bare active-space result has converged with respect
  to this sequence of orbitals" (l. 457–458) is therefore an overclaim. The open TODO at l. 468–469 already
  suspects this.

**Fix.**

- (a) Report, for every multireference run, the composition of the active orbitals (π/π*, n, σ, iodine
  character) and the SA-CASSCF natural occupation numbers. Put them in a SI table.
- (b) Re-select the active space by a reproducible rule: AVAS on the C/N/B 2p_z π system, or explicit π
  selection. Rerun at least BDP and I2-BDP at (12e,11o) with that rule. If the results move, rerun all nine.
- (c) In the sweep, report natural occupations at each size. Replace "has converged" by "does not change when
  orbitals X and Y (occupations …) are added", or remove the (14,14) row until a (16,16) step with π orbitals is
  available.
- (d) Reconcile the RHF-orbital and SA-orbital NEVPT2 numbers explicitly in the text. Until (a)–(c) are done,
  the sentence "it is not an active-space truncation error but a limitation of the second-order treatment"
  (l. 464–466) should be weakened to "is not reduced by enlarging the active space along this orbital
  sequence".

### M2. I2-BDP-0: the data argue against a near-resonant T2, but the abstract and discussion imply one

**Where:** abstract l. 52–54; sec:baseline "The distyryl failure" (l. 360–375); sec:budget l. 391–394;
sec:cas l. 452–453; Discussion l. 674; Conclusions l. 714.

**Problem.** The positions of T2 above S1 for I2-BDP-0 are:

| Method | T2 above S1 |
|---|---|
| ΔSCF-corrected pipeline | 0.51 eV |
| PBE0 | 0.44 eV |
| SA-CASSCF/SC-NEVPT2 | 0.54 eV |
| uncorrected B3LYP TDDFT | 0.08 eV |

- **Only one method dissents.** Three of the four methods agree within 0.1 eV. Only uncorrected TDDFT, which the
  paper itself shows is 0.35–0.76 eV off for S1, puts T2 near resonance. "The methods we used disagree on that
  position by about 0.4 eV" (l. 371–372) and "the position of an upper triplet that decides this rate differs
  by 0.4 eV between methods" (abstract l. 53–54) therefore misrepresent the evidence. "Decides this rate" is
  not established: the measured fast ISC remains unexplained, and the near-resonance mechanism is the
  hypothesis that the better methods do not support.
- **One dye, not a class.** "π-extended dyes, whose ISC rate it underestimates 26-fold" (abstract l. 53) uses
  the plural, but there is one measured rate. The Discussion (l. 674) generalises to "π-extended ones or
  near-resonant upper triplets".

**Fix.**

- **Abstract.** Use something like: "fails for the one π-extended dye with a measured rate, whose ISC rate it
  underestimates 26-fold. Only uncorrected TDDFT places an upper triplet close enough to S1 to explain the
  measured rate, and three higher-level treatments do not, so the mechanism remains open."
- **Main text.** Keep the sensitivity statement (Φ_T moves between 0.04 and 0.67 within half an eV) but drop
  "the methods disagree by 0.4 eV". Also list the remaining candidate explanations, each untested:
  - S1 relaxation on the distyryl torsions (adiabatic vs vertical ordering);
  - vibronic (spin–vibronic) coupling to T2;
  - a T1 channel with a larger SOC at the S1 geometry;
  - a solvent effect.
  Note that I2-BDP-0 is also the dye where Table `tab:abs` misses S1 by −0.23 eV, so the pipeline's S1 is not
  trustworthy there.
- **Discussion and Conclusions.** Replace "not for π-extended ones or near-resonant upper triplets" with "not
  for the π-extended dye tested".

### M3. The abstract conflates the class margin with the rate tolerance

**Where:** abstract l. 48–50 ("the triplet-yield class tolerates about 0.4 eV unless an upper triplet is near
resonance, in which case the tolerance falls to 0.016 eV, as for the top candidate").

**Problem.** For azaBDP_Br2 the 0.016 eV is the factor-2 *k_ISC* tolerance. Its *class* margin is 0.08 eV
(l. 402–403; `accuracy_budget.md`). The class margins of 0.01–0.02 eV belong to Br2-BDP and azaBDP_IPh2. Those
two are tight because they sit at Φ_T ≈ 0.5, not because of a near-resonant triplet. The sentence as written
gives the wrong number for the wrong reason. The Conclusions (l. 712–713) phrase it correctly as a gap
tolerance.

**Fix.** "…the triplet-yield class tolerates about 0.4 eV for dyes far from Φ_T = 0.5, but only 0.01–0.02 eV
for dyes near it. When the rate runs through a near-resonant upper triplet, as for the top candidate, the rate
tolerance falls to 0.016 eV and the class margin to 0.08 eV."

### M4. The vertical-versus-measured comparison is applied to one method only

**Where:** l. 295–297, l. 326–327 (`tab:isc` ΔE_ST), l. 440–443, l. 464–466, `tab:sanevpt2`.

**Problem.** All computed energies are vertical at the GFN2-xTB S0 geometry. They are compared with absorption
band maxima for S1 and with measured (0–0-type) singlet–triplet gaps from Ly2021. The paper invokes this
mismatch only to excuse NEVPT2 ("which may account for part of it", l. 465–466). The pipeline benefits from the
same comparison:
- its −0.055 eV absorption error;
- its 0.72 vs 0.77 eV and 0.61 vs 0.59 eV ΔE_ST.

For BODIPYs the vertical energy typically lies above the band maximum by about 0.1–0.2 eV, and the vertical
S1–T1 gap exceeds the adiabatic one because T1 relaxes more than S1. The pipeline's near-perfect agreement may
therefore partly reflect error cancellation, and the NEVPT2 errors may be smaller than reported. Either way,
the budget comparison in sec:cas (errors of "three to four times the budget") rests on it.

**Fix.**
- Estimate the vertical-to-band-maximum and vertical-to-adiabatic corrections for two dyes (BDP and I2-BDP),
  for example from TDDFT S1 and ΔSCF T1 geometry optimisations or a vibronic band-shape calculation.
- Apply the corrections to *all* methods in `tab:abs`, `tab:isc`, `tab:cas`, `tab:sanevpt2` and `tab:sweep`, or
  at least state them as an uncertainty band.
- State in the pipeline section that its agreement is with a vertical-to-measured comparison.

### M5. The decision-cost inputs mix levels of theory, and the recommended target is not costed

**Where:** sec:decision-method (l. 240–264); sec:decision (l. 591–648); `tab:decision`; abstract l. 59–66;
thesis sentence l. 64–66, Discussion l. 672–681 and Conclusions l. 719–721.

**Problem.**
- **Three levels of theory.** In p = k_ISC/K′, k_ISC is the TDDFT/ΔSCF pipeline rate, ‖χ‖ is from TDDFT, and
  α_SO is from the CAS(12,11) SOC matrix (`decision_cost.md`). The cost is for a THC encoding of the *bare*
  (12e,11o) Hamiltonian. The paper's own sec:filter shows that this Hamiltonian gives k_ISC values wrong by
  orders of magnitude. So the p the machine would actually see is not the p in the table. The authors know
  this, since Limitations l. 702 says "takes its electronic inputs from TDDFT", but the abstract presents the
  costs without the caveat.
- **Untested assumption.** The S1-preparation overlap γ² = 0.5 enters every call and is not tested.
- **Thesis not costed.** The thesis sentence recommends "excited-state energies with dynamic correlation at the
  budget accuracy, combined with a gap-resolved observable and amplitude estimation" as the useful target. No
  cost or feasibility estimate for that target is given; it is the open TODO "Option B" at l. 681. As written,
  the recommendation is an argument, not a result.

**Fix.**
- (a) Say in the abstract that the costs assume the pipeline's rates as the target values.
- (b) Give the sensitivity explicitly. Decision cost scales as p^(−1/2) and as 1/γ (amplitude amplification),
  so a factor 10 error in k_ISC changes the D1 cost by about 3.
- (c) Compute ‖χ‖ in the same CAS(12,11) space as α_SO, so that the ratio α_SO/‖χ‖, which the Discussion calls
  "the cost lever", is computed at one level of theory.
- (d) Either complete the Option B estimate or soften the thesis sentence and the Conclusions to "we argue
  that…" and move the target statement to the Discussion.

### M6. Headline ranges in the abstract, decision section and conclusions do not match `tab:decision`

**Where:** abstract l. 61–64; l. 594–597; l. 607–625; l. 686; Conclusions l. 717–719.
Source: `quantum/results/decision_cost.md`.

| Text (line) | Text says | Source gives | Fix |
|---|---|---|---|
| Abstract l. 61; l. 597 | p is 10⁻⁸–10⁻⁵ | 6.2×10⁻⁹–4.4×10⁻⁵ (ISC-active) | "10⁻⁸–10⁻⁴" or "6×10⁻⁹–4×10⁻⁵" |
| Abstract l. 62; Concl. l. 718 | a decision takes 10⁴–10⁶ calls | MLAE calls: D2 4.3×10²–1.8×10⁶; D1 3.4×10⁴–2.9×10⁶ | "10³–10⁶ calls" (or state D1 and D2 separately) |
| Abstract l. 62 | 10¹¹ Toffolis "for clear negatives" | true for light-atom negatives only (5.2×10¹⁰–5.4×10¹¹); heavy-atom negatives I-BDP 1.7×10¹³, I-BDP-0 6.7×10¹², I2-BDP-0 1.2×10¹³ | "10¹¹ for light-atom negatives, 10¹²–10¹³ for heavy-atom dyes far from the threshold" |
| l. 617 | light-atom negatives 5×10¹⁰–4×10¹¹ | BODIPY core is 5.4×10¹¹ | "5×10¹⁰–5×10¹¹" |
| l. 594, l. 686 | α_SO/‖χ‖ = 7–150 for heavy-atom dyes | I1-BDP 5.2, I3-BDP 8.6, azaBDP_Br2 6.6, azaBDP_IPh2 156; I-BDP 728 and I-BDP-0 668 if included | "5–160 for the ISC-active heavy-atom dyes" (and say I-BDP and I-BDP-0 are excluded) |
| l. 618–619 | family-wise cost grows by ln(40N)/ln 40, "a factor of 1.6 for N = 100" | ln 4000 / ln 40 = 2.25 | "2.2". The same error is in `decision_cost.py` l. 279 |
| l. 621 | decision costs are 10²–10⁴ times the 3.8×10¹⁰ of Zhou2026PRA | D1: 1.4×10²–1.6×10⁴ ✓; D2: 1.4 (azaBDP_Ph4) to 10⁴ | "D1 costs 10²–10⁴ times…; D2 for light-atom negatives is comparable to the proxy" |
| l. 624–625 | "Our totals fall in the range reported … by Eklund et al." | Eklund2026 quote 10¹³–10¹⁶ Toffolis; ours span 5×10¹⁰–6×10¹⁴ | "overlap the lower end of the range" |
| Abstract l. 60; Concl. l. 717 | 131–193 logical qubits | `thc_pennylane_resources.md`: 131–374 (rank 22), 136–379 (rank 44); costs use 158. 193 is rank 44 with optimal swap depth | "131–374 logical qubits depending on angle batching (158 in our cost model)" |

### M7. The classical validation is thin and its key ingredient, the SOC operator, is not validated

**Where:** Methods l. 174–178 (SOC) and l. 189–191 (k_ic); `tab:isc`; "Choice of F" l. 207–211; Limitations l. 698–702.

**Problem.**
- **Two dyes carry the factor-2 claim.** The claim "ISC rates within a factor of two" (abstract l. 52) and the
  choice F = 2 rest on two dyes with nearly identical chemistry (I2-BDP, I3-BDP; ratios 2.0 and 1.8), both
  overestimated. The third measured rate (I2-BDP-0) fails by 26×. No uncertainty is attached to any rate.
- **Ad hoc SOC model.** The SOC is a one-electron operator whose per-element scale factor reproduces the atomic
  ζ "in the same basis". For iodine this is applied on top of an ECP, which has no spin–orbit part unless an
  SO-ECP is used. This is not a standard method. Its only test is formaldehyde (El-Sayed rules), which has no
  heavy atom. V1 = 13.2 cm⁻¹ for I2-BDP determines both the validation and the cost tables (‖χ‖ = 1120 cm⁻¹).
- **One k_ic for all dyes.** k_ic is fixed at the BDP value for all dyes, including the aza-BODIPYs, which are
  known to have much faster non-radiative decay. This directly shifts the D2 threshold and the class margins.

**Fix.**
- Benchmark V1 (and V2 for I2-BDP-0) against a standard SOC treatment, such as SOMF/Breit–Pauli or X2C with an
  all-electron or SO-ECP basis (ORCA or PySCF), for I2-BDP and I2-BDP-0.
- State the factor-2 agreement as "for the two iodinated compact dyes with measured rates" in the abstract.
- Give the k_ISC range from the λ_s/S variation (l. 333–336) as an error bar in `tab:isc`.
- Use measured k_ic for the aza-BODIPYs, or show that D2 does not change over a plausible k_ic range.

### M8. Hamiltonian-approximation errors are absent from the budget

**Where:** sec:resources (l. 576–589); sec:decision-method l. 249–264.

**Problem.** The budget allows about 2 mHa per state (l. 387). The cost model charges 2 mHa to the single QPE
that centres the filter (l. 257–258). Two other errors are not quantified:
- the THC factorization error at rank 22 (and 44) on the S1–T_n gaps;
- the effect of the 20/27-bit precision rules of Zhou2026PRA on the filter.

Both errors enter the gap the decision depends on.

**Fix.** Report the S1–T1 and S1–T2 gaps of the BODIPY core from the THC-factorized Hamiltonian at rank 22 and
44 against exact CASCI. State what fraction of the 0.07 eV budget they use. If rank 22 is not enough, rerun the
resource estimate at the rank that is.

### M9. Novelty and related-work positioning

**Where:** Introduction l. 88–156.

Positioning is now mostly sound: the paper cites Zak2026, Nykanen2024, Genin2025, Eklund2026, Roggero2020 and
Rall2020, and l. 120 says "We do not propose a new quantum algorithm". Three points remain:
- **Last sentence of the related work.** "What has not been done for photosensitizers is … to cost the decision
  rather than the shot" (l. 153–156) holds only because of the "for photosensitizers" qualifier. Eklund2026
  costs yields and rates including amplitude estimation. Say this in the same sentence, for example: "as Eklund
  et al. did for photochemical yields".
- **Proxy comparison.** The comparison with the gap-free proxy at l. 621–624 is fair as written, but it should
  appear in the abstract too. A reader of "10¹¹–10¹⁴ Toffolis" will compare it with 3.8×10¹⁰ without the
  caveat.
- **"verified … to 0.03%".** This is a statement about QSVT phase-finding on a 6–11 wire simulation, not about
  the observable's accuracy. Add "(2% for the smallest rate)" in the abstract (l. 57). The I2-BDP CAS(4,4)
  entry of `tab:filter` is off by 2.4%.

---

## Minor points

1. **l. 601, Grover-law deviation.** The text says it is followed "to 4×10⁻⁹", but `ae_demo.json`
   (`grover_law_max_dev`) gives 5.4×10⁻¹⁰. Use the generated value.
2. **l. 596–606, where the statistics come from.** The MLAE and canonical-QAE error statistics come from the
   analytic sin² model with an ensemble over nearby amplitudes. The circuit was used to verify the Grover law
   and a few canonical points (`circuit_qpe`). Say so. "verified on the circuit" covers only the law.
3. **MLAE error floor.** At 2.7–3.3% the MLAE error stops decreasing (8×10⁵–3.3×10⁶ calls; c_eff rises from 45
   to 255). Limitations mentions 5%. Show this in Figure `fig:ae` or in the SI, and say what causes the floor,
   such as the grid of Grover depths or the shot count per depth.
4. **`ae_demo.json`, key names.** The keys `q90_*` and `c_mlae_q90` hold 95th-percentile values (CONF = 0.95
   in `ae_demo.py`). Rename them to `q95_*` before the data are archived.
5. **l. 55 and l. 432, "20–25 times the budget".** 1.4–2.0 eV over 0.07 eV is 20–29 times. Use "20–30 times" or
   divide by the 0.08 eV upper end and say so.
6. **l. 55, l. 431–432, `tab:cas`, orbital choice.** The 1.4–2.0 eV CASCI error is on RHF orbitals. With
   SA-CASSCF orbitals the same active space misses S1 by 0.65–1.05 eV for the compact dyes (CASSCF S1 in
   `references/results/*.json`). State the orbital choice in the abstract and in the `tab:cas` caption, and
   report both, because the "20–25 times" headline depends on it.
7. **l. 444, "within 0.04–0.14 eV of the pipeline".** The unrounded values give 0.045–0.140 eV, and from the
   table as printed 0.05–0.14. Use "0.05–0.14".
8. **SA weighting.** SA-CASSCF is averaged over three singlets with equal weights, and triplets are computed by
   CASCI on these singlet orbitals. T2 is then evaluated with orbitals not optimised for any triplet. Mention
   this in Methods, and for one dye compare with SA over singlets and triplets.
9. **Basis and geometry.** NEVPT2 in def2-SVP and geometries from GFN2-xTB are both light choices for 0.07 eV
   claims. Add one def2-TZVP NEVPT2 point (BDP) and one B3LYP-optimised geometry (BDP, I2-BDP) to show the size
   of these effects.
10. **l. 442–443.** "with the same sign for every dye" is correct, but it is the systematic sign the vertical
    comparison (M4) would produce. Say so.
11. **l. 364–369.** The PBE0 numbers (compact dyes within 0.00 to +0.05 eV; k_ISC 1.6×10⁹; I2-BDP-0
    7.6×10⁶, T2 0.44 eV above S1) and the conformer statement (factor below 1.5) have no generated source among
    the files checked. Add them to the SI with their output files.
12. **Methods l. 176, SOC prose.** "e.g. 28 cm⁻¹ for C, 5069 cm⁻¹ for I". Give the full table of scale
    factors and ζ values in the SI.
13. **sec:decision-method l. 236–243.** Define ρ_max, and say that K′ uses the same ρ_max as the QSVT
    normalisation of eq. `eq:filter`.
14. **sec:decision-method l. 254–257.** "a sum of 10⁴ Slater determinants" needs a justification or a citation.
    Is it the number of determinants above a weight threshold in the CAS(12,11) S1 vector?
15. **`tab:decision` caption.** Add that D1 is listed only for dyes with k_ISC ≥ 10⁶ s⁻¹, as in
    `decision_cost.md`. Also add the per-call cost column (1.6–2.1×10⁸), so readers can check "total/per call".
16. **`tab:filter`.** Give V1 with a unit-consistent precision. "0.000" and "0.001 cm⁻¹" next to "11.49" reads
    as a rounding artefact. Add a column with the relative error.
17. **Open TODOs.**
    - TOC graphic;
    - re-checking the Ly2021 values and the 5 nm aza-BDP-Ph4 discrepancy (l. 282–284);
    - Option B (l. 681);
    - acknowledgement;
    - SI contents (l. 730–731);
    - Zenodo DOI (l. 735).

    The Ly2021 check matters most, because three of the four measured gaps and all three measured rates come from
    that table.
18. **Error bars generally.**
    - Give the absorption mean error with its standard deviation and n = 6 (done).
    - Give the λ_s range for k_ISC (M7).
    - Give the cost tables with a factor range from the p uncertainty (M5).

    A referee will ask for all three.
19. **Reproducibility.** The Methods should name, for each table, the script and the command that produces it,
    as the Data Availability statement implies:
    - `run_sa_nevpt2.py` and `run_cas_sweep.py` with arguments;
    - `decision_cost.py`;
    - `ae_demo.py`;
    - `accuracy_budget.py`;
    - `tolerance_curve.py`.

    Tag the commit used for the paper.

---

## Numbers checked and found consistent

- **Absorption.** `tab:abs` mean error −0.055 eV, SD 0.021 eV. The compact-dye errors and the distyryl errors
  of −0.21 to −0.23 eV match.
- **ISC.** The 26-fold underestimate for I2-BDP-0 (4.5×10⁸ / 1.7×10⁷) is correct.
- **Accuracy budget (`accuracy_budget.md`, `tolerance_curve.json`).**
  - Factor-2: 0.062–0.083 eV.
  - Factor-1.5: 0.035–0.056 eV.
  - Factor-3: 0.12–0.14 eV.
  - Factor-10: 0.18–0.38 eV.
  - λ_s ranges: 0.03–0.07 eV and 0.09–0.13 eV.
  - Class margins: 0.4 eV for two dyes; 0.08 eV for azaBDP_Br2; 0.01–0.02 eV for Br2-BDP and azaBDP_IPh2.
  - azaBDP_Br2 k-tolerance: 0.016 eV.
- **`tab:cas`.** CASCI error 1.44–2.01 eV. NEVPT2 error 0.06–0.27 eV. The two-larger, two-comparable statement
  holds.
- **`tab:sanevpt2` (`references/results/*.sa_nevpt2_12_11.json`).** All nine dyes, all SA-CASSCF runs
  converged.
  - S1 error 0.51–0.79 eV for the compact dyes and 0.38–0.52 eV for the distyryl dyes.
  - S1–T1 0.80–1.03 eV, errors 0.22–0.27 eV.
  - T2–S1: 0.37–0.45 eV for the iodinated and phenyl dyes; 0.17 eV for aza-BDP-Ph4; 0.50–0.73 eV for the
    distyryl dyes.
- **`tab:sweep` (`sweep_partial_2026-11-17.md`).**
  - All entries match, including (14,14).
  - S1 falls 0.74 eV (BDP) and 0.22 eV (I2-BDP).
  - At (12,12): S1 0.60–0.63 eV above experiment; gaps 0.26 and 0.48 eV too large.
  - NEVPT2 gaps 1.00/0.99/0.96 and 0.82/0.87/0.88 eV.
  - T2–S1 range 0.24 eV for I2-BDP.
- **sec:resources (`thc_pennylane_resources.md`, `decision_cost.md`).**
  - 2.8–3.5×10³ Toffolis per step; 374/239/158/131 wires; 10–30% overhead; rank-44 136–379 wires.
  - About 3×10⁴ walk steps and about 10⁸ Toffolis for the filter; 1.4–2.1×10⁸ per call.
- **sec:decision.**
  - α_SO 850–12,400 cm⁻¹.
  - D1 MLAE 5×10¹²–6×10¹⁴ Toffolis, 3×10⁴–3×10⁶ calls per call cost, sampling/MLAE ratio 2×10²–2×10⁴.
  - D2 heavy-atom range 6×10¹¹–2×10¹³; Br2-BDP and azaBDP_IPh2 2.3 and 3.7×10¹⁴, at 0.94 and 1.2 times the
    threshold.
  - MLAE beats SPRT for every dye.
  - All `tab:decision` entries match.
- **AE demo (`ae_demo.json`).**
  - MLAE 18% at 1.03×10⁵ calls and 9% at 2.06×10⁵.
  - Canonical 18% at 2.6×10⁵.
  - Sampling 250% and 160%.
  - π/√p = 1.33×10³.
  - c = 45 (MLAE) and 109 (canonical).
- **Filter degree and width.** The degree fit (0.554 ln(1/ε) − 0.106) and Δ = √(2λ_s k_BT) = 0.072 eV match
  `degree_check.json`.

## Sources checked

The sources below are under `paper/`:
- `analysis/accuracy_budget.md`
- `analysis/tolerance_curve.json`
- `quantum/results/decision_cost.md`
- `quantum/results/alpha_so.json`
- `quantum/results/ae_demo.json`
- `quantum/results/degree_check.json`
- `quantum/results/thc_pennylane_resources.md`
- `quantum/results/*.filtered_isc.json`
- `references/results/*.sa_nevpt2_12_11.json`
- `references/sweep_partial_2026-11-17.md`
- `references/run_sa_nevpt2.py`
- `references/run_cas_sweep.py`
- `quantum/decision_cost.py`
- `quantum/ae_demo.py`
