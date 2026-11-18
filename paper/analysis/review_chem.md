# Referee report (electronic structure / ISC rate theory): "How Accurate Must Excited States Be to Screen Triplet Photosensitizers?" (JCTC draft)

Manuscript: `paper/main.tex` (726 lines) and `paper/si.tex` as found on 2026-11-18. Line numbers refer to
`main.tex` unless marked SI. The manuscript builds without errors or undefined references
(`latexmk`, outdir `/tmp/rev-chem`). This report covers the classical pipeline and its validation, the accuracy
budget, the I2-BDP-0 discussion, the active-space analysis and the conclusions. It does not review the QSVT and
decision-cost parts in depth. A previous report (`paper/analysis/manuscript_review.md`) is not repeated where the
current text has fixed it. Points from that report that are still open are listed in a separate section.

## Summary and recommendation

The paper asks how accurate S1–Tn gaps must be for an ISC-based screening decision on BODIPY dyes. It answers
about 0.06–0.08 eV for a factor-2 rate and 0.016 eV when a near-resonant upper triplet carries the rate. It then
compares this budget with a TDDFT/ΔSCF pipeline and with CASCI/NEVPT2 in two (12e,11o) active spaces. Since the
previous round, the paper has improved a lot:
- an AVAS-type π selection (`run_pi_cas.py`) with orbital compositions;
- an honest discussion of the canonical window;
- caveats on SOC, k_ic and the vertical comparison;
- consistent numbers.

Almost every number I checked matches its generated source (see the last section).

My recommendation is **major revision**. The budget idea is sound and useful. Three of the paper's central chemical
claims, however, rest on a more selective or less controlled reading of the data than the text admits:
1. The I2-BDP-0 paragraph now claims a consensus that T2 lies about 0.5 eV above S1. The repository's own data do
   not support that consensus.
2. "The pipeline meets the budget for compact dyes" holds only for B3LYP, and only through an uncontrolled choice
   of how the upper triplets are treated.
3. The π-only (12e,11o) NEVPT2 result is read as a failure of second-order dynamic correlation. It has a known,
   specific cause in the reference space and protocol that the paper does not test.

There is also an error in the natural-occupation logic of `run_pi_cas.py`, and the rate model has several
uncaveated or partly caveated inputs.

---

## Major points

### M1. I2-BDP-0: the "three of four methods agree" statement is not supported by the data (l. 344–362)

**Problem.** Lines 352–355 say "Three of the four methods therefore agree within 0.1 eV that T2 lies about half an
eV above S1; only uncorrected TDDFT … places it near resonance." Four problems:

- **One functional is left out.** The functional scan documented in `docs/VALIDATION.md` (l. 194–216; data in
  `validation/runs_xc/` on branch `results-tars`) includes BHandHLYP, and BHandHLYP places T2 only **+0.12 eV**
  above S1 for I2-BDP-0. Its dominant channel is S1→T2 (`runs_xc/bhandhlyp/I2-BDP-0.json`: T2 dE = +0.120 eV,
  V = 2.2 cm⁻¹). The paper reports PBE0 (l. 350–352) but not BHandHLYP. So two of the five treatments in the
  repository are near resonance, not one of four.
- **Two of the "three" are not independent.** The B3LYP and PBE0 pipelines use the same construction:
  - S1 is shifted down by ΔSCF (−0.46 and −0.42 eV);
  - all triplets, including T2, get only the T1 shift (−0.02 eV);
  - this follows from `corrections` in `validation/runs/set3/I2-BDP-0_c0.json` and the PBE0 file.

  Their agreement on T2–S1 therefore reflects the construction, not two independent determinations.
- **The third (canonical-window NEVPT2, SI Table S2) is uncontrolled.** By the paper's own Table 3, the canonical
  window of a distyryl dye holds only about 3 of the 11 dipyrrin π orbitals (BDP-0). NEVPT2 in that window misses
  S1 of I2-BDP-0 by +0.46 eV.
- **The pipeline's S1 is too low for this dye.** It is 0.23 eV below the band maximum (Table 1), and by l. 284–286
  vertical energies lie 0.1–0.2 eV above band maxima. If the pipeline's T2 is kept and S1 is placed where
  experiment puts it, T2–S1 becomes about 0.51 − (0.23 + 0.1…0.2) ≈ **0.08–0.18 eV**. That is in the range where
  the paper itself shows Φ_T moving from 0.04 to 0.67. With the pipeline's own error estimate, the evidence is
  therefore ambiguous and, if anything, leans toward near resonance.

The paper also leaves out a test that supports its framing:
- **Uncorrected TDDFT fails the compact dyes.** Raw B3LYP TDDFT gives k_ISC = 3.2×10¹¹ s⁻¹ for I2-BDP (90× too
  fast, via T2 at −0.16 eV) and Φ_T = 0.47 for BDP, where no ISC is observed (`raw_tddft` blocks in
  `validation/runs/BDP.json` and `I2-BDP.json`).
- **The shift test is unreported.** In the "higher-triplet shift test" (`docs/VALIDATION.md` l. 160–161),
  shifting T2+ like S1 breaks the compact dyes: Φ_T = 0.87 for the BODIPY core and 0.56 for BDP.

These results argue against trusting raw TDDFT for I2-BDP-0. They do not show where T2 really is.

**Fix.**
- Replace l. 352–355 with a statement of the actual spread: T2–S1 for I2-BDP-0 is 0.08 eV (raw B3LYP), 0.12 eV
  (BHandHLYP/ΔSCF), 0.44 eV (PBE0/ΔSCF), 0.51 eV (B3LYP/ΔSCF) and 0.54 eV (canonical-window NEVPT2). Say that the
  ΔSCF entries share one construction and that anchoring S1 to experiment gives about 0.1–0.2 eV. Conclude that
  the position of T2 is uncertain by about 0.4 eV, which is roughly 25 times the 0.016 eV near-resonance budget.
- Report the raw-TDDFT k_ISC and Φ_T for all dyes, for example as an extra column in Table 2, and the shift test.
  They are the real evidence that the S1-only correction is right for compact dyes.
- Determine T2 for I2-BDP-0 (and for azaBDP_Br2, l. 393–394) with one consistent method for S1 and T2. Options:
  - T2 = E_ΔSCF(T1) + [E(T2) − E(T1)] from TDDFT (an internally consistent triplet manifold);
  - ROKS/ΔSCF for T2;
  - ADC(2) or STEOM-CCSD;
  - the π-AVAS SA-NEVPT2 that is still running, on a π space that covers the styryl arms (see M5d).

### M2. "The pipeline meets the budget for compact dyes" is specific to one functional and one calibration choice (l. 279–321, 640–642, 685–687; Table 2)

**Problem.** For I2-BDP, the one dye pair that tests the factor-2 claim, the ΔSCF S1–T1 gap and the rate depend
strongly on the functional (`runs_xc`, `docs/VALIDATION.md` l. 194–216):

| | B3LYP | PBE0 | BHandHLYP | exp |
|---|---|---|---|---|
| S1–T1 (eV) | 0.61 | 0.71 | 1.14 | 0.59 |
| k_ISC (s⁻¹) | 6.9×10⁹ | 1.6×10⁹ | 1.8×10⁶ | 3.5×10⁹ |

- **PBE0 is already outside the budget.** PBE0 is the better functional for compact-dye absorption (l. 350–351).
  It misses the I2-BDP gap by 0.12 eV, beyond the 0.08 eV budget, and the rate by 2.2× (2.4× for I3-BDP), outside
  F = 2. The text quotes the PBE0 rate (l. 351) without saying either.
- **The spread is far larger than the budget.** Across three global hybrids the gap spans 0.5 eV, about 7 times
  the budget. The B3LYP agreement on two nearly identical dyes may be fortuitous.
- **The I2-BDP agreement depends on how the upper triplets are treated.** It requires T2 (29 cm⁻¹ of SOC) to be
  pushed 0.32 eV above S1. That happens only because S1 gets the ΔSCF shift (−0.50 eV) and T2 gets the T1 shift
  (−0.03 eV). The paper flags this inconsistency for azaBDP_Br2 (l. 393–394) but not for the validation dye.
- **The MLJ parameters are set a priori.** λs = 0.10, ħω = 0.17 eV and S = 0.8 (l. 164) have no citation and were
  not computed. Their stated range moves k_ISC(I2-BDP) by more than a factor of 10 (l. 319–321). A factor-2
  agreement is therefore not a discriminating test of the gaps. Also, the quoted range 2.3×10⁹–2.8×10¹⁰ comes from
  an older run (`docs/VALIDATION.md` l. 77–85, where the default row gives 8.0×10⁹, not the current 6.9×10⁹).
  Recomputed from the current `validation/runs/I2-BDP.json` over λs ∈ [0.05, 0.2] and S ∈ [0.5, 1.2] it is
  **1.7×10⁹–2.5×10¹⁰** s⁻¹. The ratio to BDP is ≥ 4×10¹¹, so "above 10¹¹" stays correct.

**Fix.**
- Qualify the claim everywhere it appears (abstract l. 47–48, l. 181–183, Discussion l. 641, Conclusions
  l. 686): "with B3LYP, for the two iodinated compact dyes with measured rates".
- Report the functional spread of ΔE_ST and k_ISC in Table 2 or the SI.
- Describe the treatment of T_n≥2 as a calibrated choice and cite the shift test.
- Justify λs, ħω and S, or compute them from S1/T1 optimizations, which l. 287 already plans.
- Update the range at l. 319–320.

### M3. The rate model and the budget leave out errors of the same size as the gap budget (l. 147–168, 170–183, 369–381)

**Problem.**
- **(a) Vertical gaps in an MLJ rate.** In the MLJ expression (eq 1) ΔE is the adiabatic (0–0) energy
  difference, with λs and S·ħω (0.24 eV in total here) carrying the reorganization. The pipeline uses vertical
  differences at the S0 geometry. The error is of the order of the difference in S1 and T1 relaxation, typically
  0.1–0.2 eV for BODIPYs (l. 285–286), so it is larger than the budget. The text discusses vertical versus
  measured only for the comparison with experiment (l. 282–287), not as an input error of the rate model.
- **(b) SOC held fixed.** The budget holds SOC fixed, but k ∝ V². A factor-2 rate tolerance is equivalent to
  ±41% in V. The SOC model is non-standard: a scaled one-electron operator on top of a scalar ECP for iodine, with
  no benchmark (l. 151–153, TODO). So the gap budget and the SOC error compete for the same factor of 2. The
  "pipeline meets the budget" claim cannot tell good gaps from compensating gap and SOC errors.
- **(c) Fixed k_r and k_ic.** The paper caveats k_ic only for the aza-BODIPYs (l. 167–168), but:
  - Measured photophysics in `validation/references.json` give an IC rate of 5.4×10⁸ s⁻¹ for I2-BDP, 3.6 times
    the assumed 1.5×10⁸.
  - The Strickler–Berg k_r is 2.0–2.9×10⁸ s⁻¹ against 1.2–1.3×10⁸ measured for BDP and I2-BDP, 1.6–2.4 times
    too high.
  - For the two near-threshold candidates (Br2-BDP at 0.94 and azaBDP_IPh2 at 1.2 times k_th = k_r + k_ic,
    l. 591–592), the class is therefore set by errors in k_r + k_ic, not by the 0.01–0.02 eV gap margins quoted
    at l. 395.

**Fix.**
- Write the budget as an error budget in ln k_ISC: gap (via ∂ln k/∂ΔE), SOC (2 ln V), and k_r + k_ic for the
  class decision. State what share of the factor 2 each takes.
- Benchmark V1 and V2 against SOMF/Breit–Pauli or X2C with an all-electron or SO-ECP basis for I2-BDP and
  I2-BDP-0 (the TODO at l. 152).
- Compute adiabatic S1 and T1 for BDP and I2-BDP and say how much the MLJ input changes.
- Redo the D2 margins of Br2-BDP and azaBDP_IPh2 over a factor-2 range of k_r + k_ic.

### M4. The budget is framed as "what the screening decision needs", but 0.07 eV is a rate tolerance (l. 640, 683–684; also l. 369–381)

**Problem.**
- **Rate tolerance versus decision tolerance.** Section 3.2 shows two different things:
  - the factor-2 *rate* tolerance is 0.06–0.08 eV;
  - the *class* decision (Φ_T above or below 0.5) tolerates more than 0.4 eV for seven of the nine dyes and
    0.01–0.02 eV only for dyes at the threshold.

  The Discussion (l. 640) and the Conclusions (l. 683) nonetheless say "the screening decision needs these gaps to
  about 0.07 eV". That merges the two.
- **Dyes with no ISC set the lower edge.** The 0.06–0.08 eV range comes from all nine validation dyes. It
  includes the BODIPY core, BDP and I-BDP, with k_ISC = 0.01–3 s⁻¹, and the range's lower edge (0.062–0.065 eV)
  comes from these dyes. For them a factor-2 accuracy in k_ISC has no decision meaning.
- **The tolerance is local.** It is computed around the pipeline's own operating point. For I2-BDP-0, where the
  pipeline's mechanism is probably wrong (M1), the relevant tolerance would be the near-resonance value
  (around 0.016 eV), not 0.079 eV.
- **The λs dependence is partly a model artefact** (l. 373–375; `accuracy_budget.md`). At λs = 0.05 eV the
  single-mode progression is resolved. The −0.61 eV T1 channel of I2-BDP and I3-BDP then sits on the flank of the
  v = 3 peak at −0.56 eV. That is why their tolerance falls from 0.083 to 0.032 eV while the dyes without iodine
  barely change. A multimode or broadened FCWD would smooth this.

**Fix.**
- State the budget per decision: D1 (rate to within a factor F) and D2 (class).
- Quote the rate tolerance over the ISC-relevant dyes (k_ISC ≳ 10⁶ s⁻¹: 0.079–0.083 eV for the validation dyes,
  0.062–0.083 eV with the candidates except azaBDP_Br2 at 0.016 eV).
- Add one sentence on the local character of the tolerance and on the vibronic structure at small λs.

### M5. The π-only (12e,11o) conclusions are well posed but over-interpreted (l. 245–256, 420–500; Tables 3 and 4; SI Tables S1–S3)

The AVAS-type selection is a real improvement. I checked `run_pi_cas.py`:
- the projector S₁₂S₂₂⁻¹S₂₁ onto MINAO 2p_z;
- the separate occupied and virtual diagonalization;
- the semicanonicalization of the core, active and virtual blocks;
- the composition analysis and the M_S = 1 CASCI/NEVPT2 with spin-free RDMs.

All of these are correct (the error in the natural-occupation logic is in M6). Four problems remain.

- **(a) Is the π-only space the right reference for a cyanine-like S1?** The BODIPY S1 is an ionic-type,
  cyanine-like ππ* state, and T1 is covalent. For ionic ππ* states, a π-only CAS reference misses the σ–π dynamic
  polarization. Second-order PT on top of it is known to overestimate such states; Angeli's analysis of the ethene
  V state is the prototype. This would produce what Table 4 shows:
  - a uniform +0.54–0.59 eV error in S1 and a uniform +0.33–0.37 eV error in S1–T1;
  - smaller gap errors in the canonical window (+0.22–0.27 eV, SI Table S2), which contains σ orbitals.

  So the result does not separate "dynamic correlation" from "reference space". The statements that the AVAS
  results "show" the canonical agreement comes from mixed character (l. 490–491) and that a larger bare active
  space is not the target (l. 646) are not established.

  *Fix.* Test the σ–π effect directly for BDP and I2-BDP. Options:
  - (12e,11o)π plus a few σ/σ* orbitals of the dipyrrin (AVAS on C/N 2s, 2p_x, 2p_y), or a RAS(SD) σ space;
  - QD-NEVPT2 or XMS-CASPT2 on the same reference;
  - an independent coupled-cluster-type reference (ADC(2), CC2/CC3 or STEOM-CCSD, as in the BODIPY benchmarks the
    paper cites: Chibani 2013, Momeni 2015, Berraud-Pache 2020).

  Rephrase l. 462–463 and l. 646 until one of these is done.
- **(b) The state-averaging protocol is not tested.** Equal-weight SA over S0–S2, with T1–T3 as CASCI roots on
  singlet-optimized orbitals (l. 249–250), puts the triplets on orbitals not optimized for any triplet. Including
  S2 also dilutes S1. The gap error this causes has unknown sign and plausibly 0.1 eV magnitude: the canonical
  and AVAS SA-NEVPT2 gaps already differ by about 0.1 eV. A second issue affects T2:
  - SC-NEVPT2 inverts T2 and T3 for the BODIPY core, BDP and I-BDP (T3 lies 0.01–0.04 eV below T2;
    `*.pi_cas_12_11.json`, `sa_casscf_orbitals.E_NEVPT2_triplets_eV`).
  - Table 4 lists the CASCI-root-2 value, +0.41 and +0.40 eV for BDP and I-BDP. The lowest upper triplet lies at
    +0.36 eV.
  - Near-degenerate roots that swap at second order call for QD-NEVPT2.

  *Fix.* For BDP, compare SA(3S) with SA over singlets and triplets (`state_average_mix`, e.g. 3S+3T) and with a
  separate triplet SA-CASSCF. Report the lowest NEVPT2 upper triplet, or QD-NEVPT2, and change l. 457 to
  "0.36–0.44 eV".
- **(c) The data base is small and partly redundant.** The headline "0.33–0.37 eV, five times the budget"
  (abstract l. 50–51, l. 454–457, l. 643–644, Conclusions l. 688–689) rests on three dyes.
  - BDP and I-BDP have the same π space and give the same result to 0.02 eV, so there are effectively two
    independent points.
  - The calculations use def2-SVP and GFN2-xTB geometries. The bond-length alternation of cyanine-like dyes
    shifts S1 by tenths of an eV.
  - Vertical energies are compared with 0–0-type gaps. On this point l. 284–286 suggests the comparison may
    inflate the errors of the other methods. For the S1–T1 gap the sign probably goes the other way: T1 of BODIPYs
    usually relaxes more than S1, so the adiabatic gap is *larger* than the vertical one, and a correction would
    *increase* the NEVPT2 gap error. (The previous report, M4, had this sign wrong.)

  *Fix.* State how Ly et al. obtained ΔE_ST (for example 77 K phosphorescence onset), compute the matching
  quantity for BDP and I2-BDP (TODO l. 287–288), and add one def2-TZVP NEVPT2 point. Until then, give the error
  as "about 0.3 eV (2 dyes, def2-SVP, vertical)" rather than "five times the budget".
- **(d) The dipyrrin space is not the chromophore π space for the distyryl dyes.**
  - The selection criterion stated at l. 426–428 (eigenvalues ≥ 0.90 kept, < 0.1 for the first rejected) already
    fails for BDP-0. The lowest kept virtual eigenvalue is 0.848, and the first rejected occupied eigenvalue is
    0.152 (`BDP-0.pi_cas_12_11.json`, `avas_occ_eigs` and `avas_vir_eigs`).
  - The π system of BDP-0 has 27 centres. The 11-centre space truncates the conjugation that makes these dyes
    different, and they are the case that matters (M1).

  *Fix.* Restrict l. 426–428 to the compact dyes and give the BDP-0 numbers. For the pending distyryl runs (TODO
  l. 464–465), either select by threshold over the full conjugated π system or state that (12e,11o) is a truncated
  space there. Otherwise those results cannot support the I2-BDP-0 discussion.

### M6. `run_pi_cas.py`: the "state-averaged natural occupations" are not natural occupations of the state-averaged density (script l. 224–242; main l. 459–461)

**Problem.** `natocc()` diagonalizes each root's 1-RDM and sorts it. `natocc_sa_singlets` is then
`np.mean(..., axis=0)` of these *sorted* vectors (l. 242). This is not the spectrum of the averaged density
γ_SA = (1/3)Σ_I γ_I, which the docstring (l. 15–17) and the text ("1.41 and 0.62 in the state average", l. 460)
claim. The roots' natural orbitals differ, so the two differ.

By Ky Fan's inequality the mean of sorted spectra majorizes the spectrum of γ_SA. The true SA occupations are
therefore *less* extreme, and the qualitative statement "no inert orbitals" survives. But:
- the quoted HOMO and LUMO values 1.41 and 0.62 are not SA natural occupations;
- even in the script's own measure the range is 0.02–0.11, not "0.02–0.09". The deviations from 2 are 0.104 for
  BDP, 0.103 for I2-BDP and 0.112 for the BODIPY core.

The inference also needs care. Occupations of 1.98–1.90 and 0.05–0.10 show weak, dynamic-like correlation of the π
orbitals. They show that no orbital is inert. They do not show that the space is an adequate reference (see M5a).

**Fix.**
- Diagonalize γ_SA. Use `mc.fcisolver.make_rdm1(mc.ci, 11, 12)` on the SA solver, or `mcscf.addons.make_rdm1s`
  averaged over the roots, or run with `mc.natorb = True` and read `mc.mo_occ`.
- Rename the current field, for example `mean_sorted_root_natocc`.
- Correct l. 460–461.

---

## Minor points

1. **Table 2, aza-BDP-Ph4.** V1 = 0.103 cm⁻¹ and k_ISC = 1.39×10⁵ s⁻¹ in `validation/runs/azaBDP_Ph4.json`. The
   table prints 0.11 and 1.5×10⁵, which are the values of an older run still in `docs/VALIDATION.md` l. 62. Use
   0.10 and 1.4×10⁵, which also matches `accuracy_budget.md`.
2. **l. 461.** "CASCI on the AVAS orbitals misses S1 by 1.2–2.3 eV" holds for the five compact dyes. For BDP-0 the
   miss is 2.66 eV (4.63 against 1.97 eV). Either name the dyes or give 1.2–2.7 eV.
3. **l. 494–497.** NEVPT2 was run only up to (12,12) (SI Table S3). "Growing the window from (8,8) to (14,14) …
   leaves the NEVPT2 S1–T1 gap flat" should read "from (8,8) to (12,12)" for the NEVPT2 part.
4. **α_SO in two "canonical windows".** Table 3 and l. 564 use RHF canonical orbitals (`run_pi_cas.py`): 8315 cm⁻¹
   for I2-BDP. Table 5 and l. 559 use B3LYP Kohn–Sham canonical orbitals (`paper/quantum/alpha_so.py`, via
   `quantum_pl.active_space`): 12006 cm⁻¹ for I2-BDP, 199 against 133 cm⁻¹ for the core. Both are called "the
   canonical (12e,11o) window". State the orbitals in the Table 5 caption and at l. 559. The factor 1.4 is itself
   evidence of how arbitrary the window's SOC content is.
5. **l. 370–371.** "Downhill side of the MLJ curve": the channels at −0.5 to −0.7 eV lie in the inverted
   (energy-gap-law) region, where k grows as the gap closes. Say so; that explains why +0.1 eV *increases* k_ISC.
6. **l. 376.** "2 mHa per state" assumes independent S1 and T_n errors. Method errors in the two states are
   usually correlated, which changes the per-state requirement in either direction. One clause suffices.
7. **l. 141, triplet manifold.** The TDA triplets are computed in a 10 eV orbital window (about 40×40 pairs;
   `meta` in the run files), while the singlets use 60 eV. Give the triplet-window convergence for T2 and T3; they
   decide the upper channels.
8. **l. 147–151, SOC with an ECP.** Say explicitly that the scalar ECP has no spin–orbit part. The iodine SOC
   comes entirely from the scaled one-electron operator acting on the pseudo-valence orbitals, whose ⟨r⁻³⟩ is not
   the all-electron value. The atomic fit corrects this only on average.
9. **Table 2 caption.** aza-BDP-Ph4 is not from Ly et al. and was not measured in toluene; `references.json` gives
   CHCl3. Its Φ_T "low" comes from secondary sources. Mark it.
10. **l. 420–431, abstract l. 49–50.** "Holds only 2.5–6.7 of the 11 dipyrrin π orbitals" is a sum of projection
    weights, not an orbital count. "Has a summed dipyrrin π weight of 2.5–6.7 out of 11" is accurate.
11. **l. 565–567.** In a strictly planar π space, the one-centre SOC between π orbitals vanishes by symmetry. The
    300 cm⁻¹ of the AVAS space comes from non-planarity (0.04 Å for I2-BDP), polarization functions and
    multicentre terms. Say this. It also implies that the Table 5 costs are for an active space that cannot carry
    the physics of the rate, so they are a lower bound for a space that includes the iodine orbitals.
12. **Table 4.** Add the CASSCF reference numbers (CASCI on SA orbitals = SA-CASSCF) and the convergence flag.
    I3-BDP has no NEVPT2 yet (TODO l. 464); the Limitations sentence "the controlled π-space results cover five
    dyes" should say "five dyes at CASSCF, four at NEVPT2, three with a measured gap".

## Points from the previous report still open

- **Vertical-to-measured corrections** (previous M4): still a TODO (l. 287–288). The caveat at l. 284–286 is new
  but gives no direction (see M5c above).
- **SOC benchmark** (previous M7): TODO at l. 152. The λs/S range is still not given as an error bar in Table 2,
  and k_ic is caveated for the aza-BODIPYs only (see M3c).
- **‖χ‖ and α_SO at one level of theory** (previous M5c): not done. It is now stated in the Table 5 caption and at
  l. 561–562, which is acceptable as a caveat.
- **THC factorization error in the budget** (previous M8): TODO at l. 552–553.
- **SA weighting** (previous minor 8): now described in Methods (l. 249–250) but not tested (see M5b).
- **Basis and geometry** (previous minor 9): mentioned only in Limitations (l. 670–671). No TZVP or DFT-geometry
  point yet.
- **PBE0 and conformer sources** (previous minor 11): the data exist (`validation/runs/set3/*_c{0,1,2}.json`,
  where k_ISC of I2-BDP-0 spans 1.2–1.8×10⁷ s⁻¹ over the conformers, and `runs_xc` on `results-tars`). They are
  not in the SI, whose contents are still a TODO (l. 708–709). Put the BHandHLYP results there too (M1, M2).
- **Re-checking the Ly 2021 values** (previous minor 17): TODO at l. 268–270. It matters most here, because all
  measured gaps and rates come from that table.

## Numbers checked

All of the following match their generated sources unless noted.

- **Table 1** (`validation/runs/*.json`, `set3/*_c0.json`):
  - λ_calc = 507/517/519/548/551/678/705/708/730 nm and the S1 energies;
  - mean error −0.055 eV for the compact dyes, distyryl errors −0.21 to −0.23 eV;
  - raw TDDFT +0.35 to +0.76 eV (core +0.76, BDP +0.63, I2-BDP +0.46, aza +0.34).
- **Table 2:**
  - ΔE_ST 0.71/0.72/0.72/0.61/0.61/0.73/0.50 eV;
  - V1 13.2 and 13.5 cm⁻¹;
  - k_ISC 2.4×10⁻², 1.1×10⁻², 2.9, 6.9×10⁹, 7.1×10⁹, 1.7×10⁷;
  - Φ_T 0.95/0.95/0.04;
  - exception: aza-BDP-Ph4 (minor 1).
- **l. 344–352:**
  - 26-fold underestimate;
  - raw TDDFT T2 at +0.079 eV with 8.0 cm⁻¹, k_ISC 1.06×10⁹ s⁻¹, Φ_T 0.67;
  - ΔSCF T2 at +0.51 eV;
  - PBE0: I2-BDP 1.64×10⁹, I2-BDP-0 7.6×10⁶, T2 +0.444 eV;
  - conformers below a factor of 1.5.
- **l. 319–321:** stale; recomputed 1.7×10⁹–2.5×10¹⁰ (M2).
- **Budget** (`accuracy_budget.md`, `tolerance_curve.json`):
  - factor 1.5: 0.035–0.056 eV; factor 2: 0.062–0.083 eV; factor 3: 0.118–0.143 eV; factor 10: 0.181–0.384 eV;
  - λs = 0.05: 0.031–0.074 eV; λs = 0.25: 0.091–0.134 eV;
  - ±0.1 eV factors 1.7–2.7;
  - class margins: more than 0.5 eV for seven dyes, 0.41–0.42 eV for aza-BDP-Ph4 and I2-BDP-0;
  - azaBDP_Br2: T2 +0.107 eV, 46.8 cm⁻¹, tolerance 0.016 eV, class margin 0.08 eV;
  - Br2-BDP and azaBDP_IPh2 margins 0.01 and 0.02 eV;
  - V1² = 1.8, 26 and 8.6 cm⁻² (l. 396–397);
  - window margins 0.15–0.21 eV (calc) and 0.00–0.07 eV (exp).
- **Table 3** (`*.pi_cas_12_11.json`):
  - Σw_π canonical 6.68/5.19/3.63/2.54/2.52/3.00;
  - n_I 0/0/2/4/5/0;
  - α_SO canonical 133/27/5579/8315/11421/11 cm⁻¹;
  - Σw_π AVAS 10.95/10.69/10.69/10.62/10.62/10.59;
  - min w_π after SA-CASSCF 0.98/0.90/0.90/0.89/0.89;
  - α_SO AVAS 8/7/80/299/304/7 cm⁻¹;
  - all SA-CASSCF runs converged;
  - the eigenvalue-gap claim fails for BDP-0 (M5d).
- **Table 4:**
  - CASCI S1 2.88/2.91/2.90/2.87/2.86; S1–T1 1.02/1.09/1.10/1.00/1.01; T2–S1 +0.23/+0.24/+0.23/+0.24/+0.23;
  - NEVPT2 S1 2.996/3.010/3.003/2.901; S1–T1 1.048/1.098/1.108/0.965;
  - errors: S1 +0.54 to +0.59 eV, S1–T1 +0.33/+0.35/+0.37 eV;
  - NEVPT2 T2–S1 0.39/0.41/0.40/0.44 by root index, 0.38/0.36/0.36/0.44 for the lowest upper triplet (M5b);
  - AVAS-orbital NEVPT2 S1 errors 0.47–0.89 eV.
- **SI Table S1** (`*.casci_ref.json`): 3.97/4.14/4.32/3.35 and 2.71/2.52/2.04/1.99 eV; CASCI errors 1.4–2.0 eV,
  NEVPT2 errors 0.06–0.27 eV.
- **SI Table S2** (`*.sa_nevpt2_12_11.json`): all 9 × 3 entries match. Also S1 errors 0.51–0.79 eV and gap errors
  0.22–0.27 eV (l. 488–489).
- **SI Table S3** (`sweep_partial_2026-11-17.md`): all entries match, including the (14,14) CASCI rows.
- **α_SO in Table 5** (`paper/quantum/results/alpha_so.json`): 39/6553/12006/12372 cm⁻¹ match, but on B3LYP
  orbitals (minor 4).
- **Measured k_r and IC rate** (`validation/references.json`): BDP 1.2×10⁸ and 1.5×10⁸; I2-BDP 1.3×10⁸ and
  5.4×10⁸, against the computed k_r of 2.85×10⁸ and 2.0×10⁸ (M3c).
