# Referee report, round 4 (chemistry and methodology)

Manuscript: "How Accurate Must Excited States Be to Screen Triplet Photosensitizers?" (JCTC), main.tex + si.tex
(si_tables.tex, si_tn_shift.tex, si_tda_full.tex, si_xc_conf.tex, si_fcwd.tex). Checked against the result JSONs in
paper/references/results, paper/analysis, paper/quantum/results, validation/runs, validation/runs_xc and
screening/runs/candidates_round1.

## Summary

The paper derives an "accuracy budget" for singlet–triplet gaps from a Marcus–Levich–Jortner (MLJ) ISC model of 15
BODIPY dyes. It tests a B3LYP/ΔSCF pipeline and (12e,11o) multireference references against measured photophysics,
and costs a fault-tolerant golden-rule filter for screening decisions. Most numbers in the abstract, main text and
tables match the JSONs; the budget, the multireference results and the decision costs are reproduced correctly. Four
problems need fixing before acceptance:
- The TDA-window check of the top candidate is internally inconsistent.
- SI Table S9 is stale.
- The SOC validation is read more favourably than the data allow.
- The fragile top candidate has not had the controls its 0.016 eV tolerance requires.

**Recommendation: major revision.** Findings 1, 2 and the text parts of 3 and 6 need no new electronic-structure
calculations. Findings 3–6 also ask for a few cheap, targeted ones.

---

## Major findings

### 1. The full-TDA check is inconsistent: the "window error" of T2 is absorbed by the ΔSCF T1 correction (major)

**Location:** main.tex:119–122, 401–404, 566–569, 696–697; si_tda_full.tex:1–21 (Table S11); script
analysis/tda_full_check.py:33–37.

**What is wrong:** The pipeline sets T_n = T_n(TDA, window) + [T1(ΔSCF) − T1(TDA, window)]. The correction for
aza-BDP-Br2 is −0.0487 eV: `corrections.T` in screening/runs/candidates_round1/azaBDP_Br2.json, which equals
E_T1 1.1230 − raw TDA T1 1.1717. The full TDA lowers T1 as well as T2. `tda_full_check.py` shifts only T_{n≥2}
(`if c["T"] >= 2`) and keeps the T1 correction computed from the windowed T1, so the T1 error is counted twice.

**Evidence** (`dE_full_minus_10_eV` in analysis/tda_full_check.json):

| Dye | ΔT1 (eV) | ΔT2 (eV) | ΔT2 − ΔT1 (eV) |
|---|---|---|---|
| aza-BDP-Br2 | −0.0201 | −0.0173 | **+0.003** |
| I2-BDP | −0.0157 | −0.0192 | −0.004 |
| BODIPY core | −0.0160 | −0.0246 | −0.009 |
| I2-BDP-0 | −0.0165 | −0.0249 | −0.008 |

Done consistently, the full TDA moves the pipeline T2–S1 gap of aza-BDP-Br2 by +0.003 eV, not −0.017 eV. Its
k_ISC then drops slightly; it does not double. The following statements are therefore artefacts:
- "the window places T2 and T3 0.016–0.025 eV too high … which doubles k_ISC of aza-BDP-Br2" (l. 120–121);
- "just above the factor-2 tolerance … doubles k_ISC (Φ_T 0.97 to 0.99)" (l. 401–403);
- the 0.017 eV window term of the near-resonance error sum (l. 567–568). The correct sum is
  √(0.003² + 0.027² + 0.003²) ≈ 0.027 eV rather than 0.032 eV. This is still above 0.016 eV, so that conclusion
  survives.
- the same statement in the Limitations (l. 696–697).

The claim that T1 "the window does not affect" (l. 122) is also literally false: TDA T1 moves by −0.016 to
−0.020 eV. It is harmless only because ΔSCF replaces T1.

**Fix (text + re-run of an existing post-processing script, no new QC):**
- Apply ΔT_n − ΔT_1 in tda_full_check.py.
- Regenerate Table S11.
- Rewrite l. 120–122, 401–404, 566–569 and 696–697. Expected wording: "the window changes T_n − T_1 spacings by
  −0.009 to +0.003 eV and no rate by more than X%".

### 2. SI Table S9 is stale and contradicts main Table 3 and decision_cost.json (major, consistency)

**Location:** si_tables.tex:71–102 versus main.tex:648–655 and quantum/results/decision_cost.json and
decision_table.tex.

**Evidence:**

| Dye | Quantity | main Table 3 = JSON | SI Table S9 |
|---|---|---|---|
| I2-BDP | D1 MLAE | 1.8e14 | 1.6e14 |
| I2-BDP | D1 canonical | 1.6e15 | 1.4e15 |
| I2-BDP | D1 channel | 1.1e13 | 1.2e13 |
| I2-BDP | D2 | **8.8e13** | **1.8e13** |
| I2-BDP-0 | D2 | 4.5e13 | 1.4e13 |
| Br2-BDP | D2 | 9.5e13 | 2.8e14 |
| aza-BDP-Br2 | D2 | 3.3e12 | 6.4e11 |
| BDP | D2 | 3.1e11 | 9.6e10 |
| BODIPY core | D2 (JSON) | 2.2e12 | 5.3e11 |

In the SI, BDP D2 by MLAE (9.6e10) is lower than the main-text "best estimator" value (3.1e11), which is impossible
under one model. The SI values are close to, but not equal to, the oracle-margin entries of the current JSON (for
example I2-BDP 1.9e13). `analysis/make_si_tables.py` reads the current JSON, so si_tables.tex was simply not
regenerated after the indifference-zone D2 model was introduced. The main-text range "D2 for light-atom dyes
1e11–2e12" (l. 603) cannot be read off the SI table as it stands.

**Fix (text-only):** Re-run make_si_tables.py and recompile. Also state in the Table S9 caption that D2 is the
indifference-zone test with no known margin.

### 3. The SOC benchmark does not support the factor-2 validation as phrased; the confound and the operator conclusion need fixing (major)

**Location:** main.tex:46–48 (abstract), 133–147, 303–305, 375–377, 706–709; si.tex:249–272 (Table S6);
references/results/I2-BDP.soc_benchmark.json and I2-BDP-0.soc_benchmark.json; psscreen/tdwin.py:229–235;
psscreen/soc.py:118–128.

**(a) The best-calibrated coupling gives about 5 times the measured rate.**
- "eff., all-electron" equals SOMF × 0.933 for every iodine-dominated channel: I2-BDP 0.9323, 0.9323, 0.932 (T5);
  I2-BDP-0 0.931, 0.931.
- 0.933 is exactly ζ_exp(I)/ζ_SOMF(I) = 5068.8/5433.3. So the 7% agreement (l. 143–145) is the atomic
  calibration ratio. It shows that a one-centre operator scaled to the atom reproduces molecular SOMF; it is not an
  independent accuracy test.
- By the authors' own reasoning, the best estimate is the recalibrated operator on all-electron states:
  V1 = 20.7 cm⁻¹ instead of 13.2 cm⁻¹. That gives k_ISC(I2-BDP) ≈ 6.9e9 × (20.7/13.2)² = 1.7e10 s⁻¹, 4.9 times the
  measured 3.5e9.
- With the parameter-free multimode line shape the same coupling gives 1.0e10 (gas phase) to 2.7e10
  (λ_out = 0.10 eV) s⁻¹.
- The "within a factor of 2" in the abstract (l. 47) and Conclusions (l. 706) therefore holds only for the ECP/RPA
  states, which the paper itself identifies as the source of the error (l. 145–146).

**(b) The attribution "comes from … the ECP and basis" (l. 145–146) is confounded.**
- The pipeline singlet is RPA: its X vector is kept, Y dropped and X renormalised (tdwin.py:229–235). The benchmark
  singlet is TDA.
- The two S1 states differ by 0.22 eV in energy: 2.77 eV RPA (pipeline raw) against 2.98 eV TDA (all-electron),
  from the JSON. Part of the V1 difference may come from RPA versus TDA rather than from the ECP or basis.

**(c) Conclusions l. 708, "the rate changes by a factor of 3–40 between spin–orbit operators", is not supported.**
- The data give 2.8 (SOMF versus effective operator, I2-BDP T1), 2.5 (recalibrated) and 5.2 (I2-BDP-0 T2).
- Nothing gives 40. The only ~40 ratio in the paper is 2.6e9/6e7 (l. 309–310), which mixes line shapes with
  operators.

**(d) Two errors in the T5 remark (l. 138).**
- "T5, 2 eV above S1" is wrong. T5 lies 1.54 eV (I2-BDP: 3.80 − 2.26) and 1.13 eV (I2-BDP-0: 2.83 − 1.70) above
  ΔSCF S1, or 0.8–1.0 eV above the raw TDDFT S1.
- "its V² differs by a factor of 1.7" holds only for I2-BDP. For I2-BDP-0 it is (29.0/3.0)² ≈ 90 (Table S6).
- Mapping channels by index between two independent calculations should be checked by state overlap, or at least
  by character.

**(e) The functional dependence also mixes gap and SOC (validation/runs_xc).** V1(I2-BDP) is 13.2 (B3LYP), 10.2
(PBE0) and 2.8 cm⁻¹ (BHandHLYP), a factor of 22 in V². The text (l. 303–304, 691) presents the PBE0 and BHandHLYP
rates as a gap effect.

**Fix:**
- *Text:*
  - Restate the I2-BDP result in the abstract and Conclusions as: factor 2 with ECP/RPA couplings, factor 2.9–5 with
    all-electron couplings.
  - Correct l. 138 and l. 708.
  - Say that the 7% agreement is the atomic calibration ratio.
  - Decompose the functional dependence into its gap and V1 parts.
- *New calculation (cheap):*
  - Evaluate the effective operator on ECP/def2-SVP **TDA** singlets, and on ECP RPA singlets with X+Y or X−Y, to
    separate the state model from the ECP/basis.
  - Run one ECP calculation with a larger iodine basis (def2-TZVP) to separate the basis.

### 4. The top candidate (aza-BDP-Br2) has not had the controls its 0.016 eV tolerance requires; "the class holds (Φ_T = 0.80–0.99)" overreaches (major)

**Location:** main.tex:393–412, 566–569, 695; Table S10, Table S11; screening/runs/candidates_round1/azaBDP_Br2.json.

**(a) An anomalous near-degenerate S2.** In the raw TDDFT, aza-BDP-Br2 has S2 only **0.10 eV** above S1 (2.245 and
2.347 eV, f = 0.66 and 0.04). In every other dye with stored singlets the S1–S2 gap is 0.36–0.52 eV (my tabulation
of `E_singlets_eV`).
- A dark state this close to S1 suggests a charge-transfer or Br-lone-pair state. Its triplet partner is a natural
  candidate for the T2 that carries V2 = 46.8 cm⁻¹. That is the largest T1–T3 coupling in the whole set, larger than
  in any iodinated dye (I2-BDP V2 = 29.4), although ζ(Br)/ζ(I) = 0.48.
- If T2 has CT or n_Br→π* character, the key assumption that T2 "errs like T1" (l. 398–399) is least justified
  exactly here. B3LYP errors for CT states are not those of a ππ* T1.
- It also raises the question of whether ΔSCF/MOM converged on the same diabatic S1 as TDDFT S1, which is the state
  used for the SOC.

**(b) SOC benchmark.** The SOC operator was benchmarked only for iodine dyes. Bromine, which carries the candidate's
rate, is unbenchmarked.

**(c) Controls tested one at a time.** Each robustness test (no shift: Φ_T = 0.80; full TDA; S1 shift) changes one
factor. Combinations were not tested. For example, the no-shift gap together with V2 three times smaller would give
k_ISC ≈ 1.6e8 s⁻¹ and Φ_T ≈ 0.3, a class flip.

**(d) Not tested at all, although each could move T2–S1 by more than 0.016 eV:**
- geometry: the DFT S0 minimum moves the BDP and I2-BDP gaps by 0.01–0.07 eV (adiabatic JSONs);
- basis (def2-SVP only);
- DFT integration grid ("level 1 for the response kernel", l. 115);
- functional (no PBE0 or BHandHLYP run for the candidates);
- T2 relaxation. The aza-BDP-Br2 T1 optimisation in references/results/azaBDP_Br2.adiabatic_T1.log is unfinished
  after 22 steps and concerns T1, not T2.

**(e) Decision-cost assumption.** The decision-cost model assumes a 0.5 eV S1–S2 gap for the S1 eigenstate filter
(SI l. 173–175). For this dye the TDDFT gap is 0.10 eV, which the paper's own sensitivity row puts at
C_S1 = 1.3e8 Toffolis.

**Fix:**
- *Text:* Report the T2 and S2 character (NTOs) and the S1/S2 near-degeneracy. Say that the class is robust only to
  single-factor perturbations. Use a per-dye S1–S2 gap in the cost model.
- *New calculation (moderate):*
  - all-electron or SOMF SOC for aza-BDP-Br2 and Br2-BDP;
  - def2-TZVP and grid-level checks of the T2–S1 gap;
  - PBE0 and ωB97X-D (or another range-separated functional) for this dye;
  - T2 relaxation from TDA gradients (run_vertical_gradients.py already exists) to get λ(T2) and the S1→T2 line
    shape;
  - ideally the single-reference check of finding 5 on this dye.

### 5. There is no reliable reference for the decisive gaps; the NEVPT2 cause analysis is incomplete (major)

**Location:** main.tex:346–358, 464–485, 673–681; Table 3; Table S3.

**What is wrong:**
- The two gaps that decide the paper's narrative are T2–S1 of I2-BDP-0 (the "distyryl failure") and T2–S1 of
  aza-BDP-Br2 (the top candidate). Each is computed only with TDDFT-derived schemes, spread over 0.08–0.72 eV
  (l. 352–354), and with SC-NEVPT2 schemes the authors themselves show to be 0.3–0.6 eV off.
- No method that is known to work for BODIPY S1/T_n spacings is applied, such as ADC(2)/SCS-ADC(2), CC2, EOM-CCSD
  or DLPNO-STEOM-CCSD. These are feasible for 40–60-atom dyes in def2-SVP or def2-TZVP. Without one, the paper
  cannot say whether the pipeline, the NEVPT2 or neither is closer for T2–S1, and "none of our methods explains"
  (l. 355) remains open-ended.
- The cause analysis (l. 482–483) lists the IPEA shift as a possible reason why CASPT2 is lower than SC-NEVPT2. The
  standard IPEA shift (0.25 Ha) *raises* CASPT2 excitation energies, so it argues the other way.
- Two more likely contributors are not tested:
  - triplets computed as CASCI roots on orbitals averaged over singlets only. The 3S+2T "samix" check moves the gap
    by only 0.02–0.03 eV, so this one is largely excluded; say so.
  - strong versus partial or quasi-degenerate contraction. QD-NEVPT2 or PC-NEVPT2 is available in PySCF or ORCA for
    BDP and I2-BDP.
- Table 3 caption: SC-NEVPT2 interchanges T2 and T3 for six dyes, not three. The interchange also occurs in BDP-0,
  I-BDP-0 and I2-BDP-0 (for example BDP-0 T2–S1 0.786 > T3–S1 0.728, from BDP-0.pi_cas_12_11.json). It does not
  occur for I2-BDP, I3-BDP or aza-BDP-Ph4. Large state-specific NEVPT2 corrections (T3 drops by about 0.4 eV relative
  to CASCI) deserve a comment.

**Fix:**
- *Text:* Correct the IPEA argument and the Table 3 caption.
- *New calculation:* An ADC(2) or STEOM/EOM-CCSD benchmark of S1, T1, T2 and T3 for BDP, I2-BDP, I2-BDP-0 and
  aza-BDP-Br2, and one QD- or PC-NEVPT2 control for BDP and I2-BDP.

### 6. The headline 0.08 eV budget is a property of the resolved single-mode progression; use the multimode line shape that is already available (major, mostly text)

**Location:** main.tex:44, 158–166, 307–319, 363–374, 704; si_fcwd.tex:105–127; analysis/tolerance_curve.json;
analysis/fcwd_multimode.json.

**(a) The single-mode tolerance is the width of one resolved line.**
- With λ_s = 0.10 eV, ħω = 0.17 eV and kT = 0.0257 eV, each vibronic line has σ = √(2λ_s kT) = 0.072 eV, while the
  lines are 0.17 eV apart. The density of states is therefore a resolved comb.
- The authors note (l. 371–373) that the I2-BDP factor-2 tolerance equals the half-width of the v = 3 line,
  0.084 eV.
- tolerance_curve.json confirms that the λ_s dependence is non-monotonic: core 0.074, 0.052, 0.069 eV and
  aza-BDP-IPh2 0.074, 0.051, 0.070 eV at λ_s = 0.05, 0.075 and 0.10 eV. A physical density of states does not do
  this.

**(b) The available multimode model gives a different gas-phase budget.** I re-evaluated their Kubo line shape for
I2-BDP on a 5 meV grid (fcwd_multimode.rho_multimode with the stored Hessian and ΔSCF minima; nothing in the
repository was written):

| Line shape | Factor-2 tolerance (eV) | k(+0.1)/k | k(−0.1)/k |
|---|---|---|---|
| Multimode, gas phase | **0.050** | 2.7 | 0.58 |
| Multimode, λ_out = 0.05 eV | 0.080 | — | — |
| Multimode, λ_out = 0.10 eV | 0.100 | — | — |
| Single-mode MLJ | 0.085 | — | — |

The budget survives with an outer-sphere term but tightens to 0.05 eV in the gas phase. Physically, λ_out for a
locally excited ππ*→ππ* transition in toluene should be small, so the gas-phase limit is not academic. Section 3.2
should report this rather than say "would be smoother and give a different value".

**(c) The bracket claim is misleading.** si_fcwd.tex:127 says the two computed reorganization energies (0.12 eV from
energies, 0.19 eV harmonic) "bracket the line-shape uncertainty". The a priori λ_s + Sħω = 0.24 eV used for every
rate and budget lies *outside* that bracket. This should be stated wherever the a priori parameters are defended
(l. 158–159). The harmonic λ(S1) and λ(T1) (0.108 and 0.127 eV) also exceed the energy-based relaxation energies
(0.078 and 0.084 eV, I2-BDP.adiabatic.json) by 40–50%. That points to rectilinear-projection or Duschinsky problems,
which deserve one sentence.

**Fix:**
- *Text plus a cheap new run of an existing script:* Add multimode tolerances (λ_out = 0, 0.05, 0.10 eV) for I2-BDP
  and BDP to Section 3.2 and Table S8. Qualify the abstract's 0.08 eV as "0.05–0.10 eV with a multimode line shape".
- Note that fcwd_multimode.json stores only the λ_out = 0 run. The λ_out rows of Table S15 are regenerated in memory
  by fcwd_table.py and are not archived.

---

## Minor findings

7. **Table S7 lists unshifted TDA T2/T3** (si_tables.tex:7–9, 29). Readers who compute T2–S1 from Table S7 get
   1.96 − 1.80 = 0.16 eV for aza-BDP-Br2 instead of the pipeline's 0.11 eV (and BDP 0.51 instead of 0.48). The
   caption says "TDA T2 and T3" but not that the T1 ΔSCF shift is excluded. *Text-only:* list the shifted values or
   say so explicitly.

8. **Sensitivity of aza-BDP-Br2 understated** (main.tex:363–364, 395). Line 363 says "a 0.1 eV error changes k_ISC
   by up to a factor of 2.7". That is true only without aza-BDP-Br2, which Table S8 gives as ×21 and ×0.013 (a
   factor of 77). Line 395, "a factor of about 20", gives only the favourable sign. *Text-only.*

9. **I2-BDP-0 alternatives use a coupling the benchmark rejects** (main.tex:347–349, 405–409). The uncorrected-TDDFT
   rate (1.1e9 s⁻¹, Φ_T = 0.67) and the "S1-shift" result (Φ_T = 0.73) both use V2 = 8.0 cm⁻¹. The all-electron
   benchmark gives 3.5 (SOMF) or 3.2 cm⁻¹ (recalibrated), which reduces both to about 2e8 s⁻¹ and Φ_T ≈ 0.3. Say
   so. It strengthens, not weakens, "none of our methods explains" this dye. *Text-only.*

10. **Like-for-like experimental gap** (Table 1, main.tex:297–299, 324). The measured ΔE_ST combines room-temperature
    fluorescence in toluene with 77 K phosphorescence in an EtI glass. The matching computed quantity (0.82 and
    0.61 eV, which I verified from the adiabatic JSONs) is in the text but not in Table 1. Add a column or footnote,
    and note the temperature and medium mismatch. *Text-only.*

11. **"At least 0.90" and Table 2** (main.tex:433 versus 454–455). The text says the AVAS orbitals have π weights of
    at least 0.90 for the compact dyes. That is correct for the AVAS orbitals (minimum 0.906 for I2-BDP and I3-BDP),
    but Table 2's column is the SA-CASSCF minimum, 0.89. Name the column consistently. *Text-only.*

12. **(14,14) step** (Table S4). CASCI(14,14) reproduces (12,12) to 0.003 eV for both dyes
    (sweep_partial_2026-11-17.md), so the added orbitals are inert and saturation is not demonstrated. The I2-BDP
    CASCI T2–S1 also oscillates (+0.16, +0.74, +0.05 eV) across sizes, a sign of root switching. Comment on this or
    drop the row. *Text-only.*

13. **Functional table** (Table S12). Add V1 per functional (B3LYP 13.2, PBE0 10.2, BHandHLYP 2.8 cm⁻¹ for I2-BDP)
    so that gap and coupling effects can be separated. See 3(e). *Text-only; data exist in validation/runs_xc.*

14. **Notation** (main.tex:187). "(2 δ ln V)²" uses δ, already defined as the gap error. Write
    (2 Δ ln V)² + (∂ ln k/∂ΔE · δ)². *Text-only.*

15. **Supporting Information list** (main.tex:725–730). It omits the multimode Tables S14–S15. The Zenodo DOI is
    still a \todo (l. 737–738). *Text-only.*

---

## Numbers verified (no action needed)

- **Absorption:** mean error −0.055 eV and SD 0.021 eV (Table S1). Uncorrected TDDFT error 0.34–0.75 eV.
- **Adiabatic gaps:** BDP 0.72 → 0.75 eV and I2-BDP 0.61 → 0.60 eV. Four-point reorganization 0.09 and 0.12 eV.
  Emission-based gaps 0.82 and 0.61 eV.
- **Rates against experiment:** PBE0 factor 2.2–2.4. SOMF 2.0e10 s⁻¹ = 5.6 times the measured rate.
- **SOC sums in Table S6:** 1135, 1108, 963, 844, 1007 and 874 cm⁻².
- **Budget ranges:**
  - factor 1.5: 0.035–0.056 eV;
  - factor 2: 0.079–0.083 eV (ISC-active dyes) and 0.062–0.083 eV (all dyes except aza-BDP-Br2);
  - factor 3: 0.12–0.14 eV;
  - factor 10: 0.18–0.38 eV;
  - over λ_s: 0.03–0.07 and 0.09–0.13 eV;
  - class margins 0.42 and 0.40 eV;
  - window margins.
- **AVAS NEVPT2 (Table 3):** errors 0.54–0.59 eV (S1) and 0.33–0.39 eV (gaps); the TZVP, DFT-S0 and samix shifts.
- **Upper-triplet treatment (Table S10):** all entries, including the factor 8.7–9 and 90 times the measured rate.
- **Main Table 3:** p = γ²(V_eff/α_SO)² and every value against decision_cost.json.
- **Decision-cost ratios:**
  - D2 MLAE/sampling 3.6–274;
  - canonical/sampling 1.2–9.4;
  - margin penalty 2–5;
  - channel-route gain 12–29;
  - 6 years for I2-BDP D1;
  - α_SO/V_eff = 150–1.3e4;
  - α_SO/‖χ‖ = 11 and ‖χ‖/V_eff = 270.
