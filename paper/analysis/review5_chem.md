# Referee report, round 5 (chemistry and methodology)

Manuscript: "How Accurate Must Excited States Be to Screen Triplet Photosensitizers?" (JCTC, revised), main.tex + si.tex
(si_tables, si_tn_shift, si_tda_full, si_xc_conf, si_fcwd). Checked against paper/analysis/*.json and *.py,
paper/references/results, validation/runs(_xc), screening/runs/candidates_round1 and paper/quantum/results. All
re-evaluations were done in /tmp/claude-0/review5c by importing the authors' own functions (fcwd_multimode.analyse,
fcwd_multimode.rho_multimode, fcwd_table.aza_t2_rows, tolerance_multimode.tolerance). No repository file other than
this report was written.

## Summary

The revision fixes the two consistency errors of round 4: the full-TDA check is now done relative to T1 (Table S11
reproduces), and Table S9 has been regenerated. It also adds a multimode Kubo line shape (SI S12) whose numbers
reproduce exactly from the scripts. The new analysis of the top candidate is the weak point. The T2 relaxation of
aza-BDP-Br2 (−0.038 eV) is a small difference of two TDA-gradient relaxation energies. One of them, S1, is twice the
ΔSCF value that the paper reports for the same dye. With the ΔSCF value the adiabatic T2–S1 gap moves from +0.07 to
−0.01 eV. In the gas-phase multimode line shape, the pipeline's own vertical gap, or the no-shift gap, puts Φ_T at
0.46–0.57. So the statement "the class holds (Φ_T = 0.80–1.00), the rate is fragile by a factor of about 50" is
not supported. The other open items are mostly carry-overs that need text changes: the near-resonance tolerance is
taken from proxy modes and a "whatever the line shape" argument that the paper's own aza-BDP-Br2 line shape
contradicts; the SOC confound is untested; and several round-4 minor points remain.

**Recommendation: major revision.** None of the new problems needs a large calculation. Most are text fixes or
re-runs of the existing scripts on files that already exist (the aza-BDP-Br2 Hessian and gradients). The exceptions
are a state-identity check (NTOs or overlaps) of the TDA S1 and T2 used for the gradients, and, as before, a
correlated single-reference benchmark of the decisive gaps.

---

## Status of the round-4 findings

| # | Round-4 finding | Status | Comment |
|---|---|---|---|
| 1 | Full-TDA check counted the T1 window error twice | **Resolved** | tda_full_check.py:51–55 now applies ΔT_n − ΔT_1. Table S11 reproduces the JSON: T2–S1 shifts from −0.010 to +0.004 eV, k(aza-BDP-Br2) −11%. Text at main.tex:122–125, 415–416 and 716 and the quadrature sum (0.027 eV) are corrected. No new error introduced. |
| 2 | SI Table S9 stale | **Resolved** | All D2 entries match decision_cost.json (for example I2-BDP 8.8e13, aza-BDP-Br2 3.3e12, BDP 3.1e11). Small caption issue: the caption says "D2 … by MLAE" without "indifference zone, no known margin" (si_tables.tex:72–74). For aza-BDP-Ph4 and BDP-0 the main-text "best estimator" is sampling (1.3e11, 1.6e11), not the listed MLAE value. |
| 3 | SOC validation over-read; RPA/TDA confound; "3–40"; T5 remark; functional V1 | **Partially resolved** | Fixed: the 7% is now called the recalibration ratio (l.148–150); the T5 energies and factors (l.141); "3–40" replaced by "about 3" (l.727–728); the RPA-versus-TDA confound is now named (l.151). Not done: the requested ECP-TDA and X±Y control calculations, so the attribution to "the states" is still untested. Abstract l.47 still says "within a factor of 2" without the all-electron counterpart (×4.9–5.6). Functional dependence is still not decomposed (see 3(e) and #13). l.390 still says "the two SOC operators … differ by 68% in V1", although l.149–151 concludes that the difference comes from the states, not the operator. |
| 4 | Controls for the fragile top candidate | **Partially resolved; new errors introduced** | T1 optimisation now converged (adiabatic.json). T2 relaxation added through vertical gradients (SI S12). Not done: S1/S2 near-degeneracy and the T2/S2 character (no NTOs; the 0.10 eV RPA S1–S2 gap is not mentioned anywhere); Br SOC benchmark (acknowledged only in the Limitations, l.717); basis, grid and functional for the candidate; combined perturbations; per-dye S1–S2 gap in the cost model (si.tex:174 still assumes 0.5 eV). The new T2 analysis brings new problems: new findings 1–3. |
| 5 | No correlated single-reference reference; NEVPT2 cause analysis | **Partially resolved** | IPEA argument removed. Table 3 caption corrected (all dyes except I2-BDP, I3-BDP, aza-BDP-Ph4). State-averaging (3S+2T) control reported (l.494–496). No ADC(2)/CC2/EOM/STEOM and no QD/PC-NEVPT2: acknowledged only as a limitation (l.716–717). The large state-specific T3 NEVPT2 corrections are not commented on. |
| 6 | Budget is a property of the single-mode comb; use the multimode line shape | **Largely resolved** | Multimode tolerances are in the text (l.384–388) and match tolerance_multimode.json: I2-BDP 0.046/0.079/0.098 eV, BDP 0.036/0.069/0.068 eV. The bracket claim is withdrawn (si_fcwd.tex:37). Remaining: the abstract's "0.05–0.10 eV" omits BDP's 0.036 eV; fcwd_multimode.json is still not an archive of the reported runs (new finding 6); the aza-BDP-Br2 tolerance still uses I2-BDP proxy modes (new finding 2). |
| 7 | Table S7 lists unshifted T2/T3 | **Resolved** | The caption now says "raw TDA, before the shift". |
| 8 | Sensitivity of aza-BDP-Br2 understated | **Partially resolved** | l.406 now gives 21 and 77. l.374–375, "a 0.1 eV error … changes k_ISC by up to a factor of 2.7", is still unqualified (aza-BDP-Br2: ×77). |
| 9 | I2-BDP-0 alternatives use the rejected V2 = 8.0 cm⁻¹ | **Not addressed** | l.358–361 (Φ_T = 0.67, "close to experiment") and l.424–425 (Φ_T = 0.73) still use V2 = 8.0. With V2 = 3.5/3.2 cm⁻¹ from SOMF/eff.-all-electron, both fall to about 2e8 s⁻¹ (Φ_T ≈ 0.3). |
| 10 | Like-for-like experimental gap | **Mostly resolved** | 0.82/0.61 eV given in the text (l.309–310). The Methods describe the 77 K EtI-glass caveat. Table 1 still compares vertical gaps only; a footnote would suffice. |
| 11 | "At least 0.90" against the Table 2 column | **Resolved** | The column is now labelled "min w_π (SA-CASSCF)". |
| 12 | (14,14) row is inert; I2-BDP CASCI T2–S1 oscillates | **Not addressed** | si.tex:138, 142: no comment on the (14,14) row or on the root switching (+0.16/+0.74/+0.05 eV). |
| 13 | V1 per functional in Table S12 | **Not addressed** | validation/runs_xc: V1(I2-BDP) = 13.24 (B3LYP), 10.23 (PBE0), 2.81 cm⁻¹ (BHandHLYP). Of PBE0's factor-4.3 drop, a factor 1.7 is coupling; for BHandHLYP the coupling alone accounts for a factor 22. |
| 14 | δ used twice in the error-budget inequality | **Not addressed** | main.tex:192 still reads (2 δ ln V)². |
| 15 | SI list; Zenodo DOI | **Partially resolved** | Tables S14–S15 are now listed (l.736). The Zenodo \todo is still there (l.745). |

---

## New findings

### 1. The T2 relaxation of aza-BDP-Br2 is not determined to better than about 0.08 eV, five times the tolerance it is meant to resolve (major)

**Location:** main.tex:302–304, 416–418, 713–714; si_fcwd.tex:39–41; fcwd_table.py:22–44;
references/results/azaBDP_Br2.vertical_gradients.json and .npz; azaBDP_Br2.adiabatic.json.

**What is wrong:** The adiabatic correction of the decisive gap, T2–S1 from 0.107 to 0.069 eV, is
λ(T2) − λ(S1) = 0.191 − 0.153 = 0.038 eV, a difference of two vertical-gradient (TDA) relaxation energies. The
gradient model is not reliable at the 0.02 eV level for this dye, and it was combined with a vertical gap from a
different geometry and a different S1 method.

**Evidence:**
- *The S1 relaxation of the gradient model is twice the ΔSCF value the paper reports for the same dye.*
  - azaBDP_Br2.adiabatic.json: E_S1(S0 min) − E_S1(S1 min) = 0.0764 eV. The main text quotes it as 0.08 eV (l.301).
  - The TDA gradient gives λ(S1) = 0.153 eV. T1 agrees (0.086 against 0.094).
  - The same bias appears in BDP (0.105 against 0.051) and I2-BDP (0.163 against 0.108), so it is systematic.
  - With the ΔSCF S1 relaxation and the TDA T2 relaxation, T2 relaxes 0.115 eV more than S1, and the adiabatic gap
    becomes **−0.008 eV**, not +0.069 eV. The paper's S1 is the ΔSCF state, so this combination is at least as
    defensible as TDA–TDA.
  - The two choices bracket the adiabatic gap at −0.01 to +0.07 eV. That spread is 5 times the 0.016 eV tolerance,
    so the gradient model cannot place T2 at the level the argument needs. Re-evaluated with the authors' code, the
    S1→T2 rate then spans 3.2e9–8.1e12 s⁻¹ (gas phase).
- *The harmonic projection breaks down for this dye.* fcwd_multimode.json (molecules[0].displacements) projects
  the optimised ΔSCF minima on the S0 modes and gives λ(S1) = 0.397 eV and λ(T1) = 0.289 eV. The energy-based
  values are 0.076 and 0.094 eV, so the harmonic model overestimates them five- and threefold.
  - In my decomposition, 0.23 and 0.08 eV of λ(S1) sit in two modes at 259 and 189 cm⁻¹. That is a large-amplitude,
    anharmonic coordinate, plausibly phenyl or BF2 motion.
  - The paper reports the 1.5× harmonic excess for BDP and I2-BDP (si_fcwd.tex:37) but not this 5× failure for the
    one dye where the multimode line shape is decisive.
- *Geometry and method mismatch.*
  - The vertical gap +0.107 eV is the pipeline value at the xTB geometry: windowed density-fitted TDA T2 shifted
    by the ΔSCF T1 correction, and ΔSCF S1.
  - The relaxation comes from gradients at the B3LYP S0 minimum. Rebuilt the pipeline way at that geometry
    (vertical_gradients.json: TDA T2 1.9462 and T1 1.1432 eV; adiabatic.json: ΔSCF T1 1.1131 and S1 1.7997 eV),
    the vertical gap is 1.9462 + (1.1131 − 1.1432) − 1.7997 = **+0.116 eV**. The consistent adiabatic estimate is
    then +0.079 eV, not +0.069 eV.
  - In the gas-phase line shape this alone halves the rate (1.6e9 against 3.2e9 s⁻¹).
- *State identity is unverified.*
  - In the raw RPA spectrum S2 lies only 0.10 eV above S1 (2.245 eV, f = 0.66; 2.347 eV, f = 0.04;
    candidates_round1/azaBDP_Br2.json). The TDA roots used for the gradients are 0.125 eV apart (2.419 and
    2.545 eV).
  - The S1 and T2 gradients are almost collinear: cos = 0.97 for the mass-weighted, projected gradients, mutual λ = 0.014 eV. The S1
    and T1 gradients are almost orthogonal (cos = 0.34; BDP −0.13; I2-BDP 0.09).
  - Near-collinear gradients are what one expects when the TDA "S1" and T2 have the same orbital configuration. By
    El-Sayed's rule that same pair would then have a small coupling, which conflicts with V2 = 46.8 cm⁻¹ computed
    between the RPA S1 and T2.
  - Nothing in the repository shows that the TDA root 1 at the DFT geometry is the bright RPA S1 used for the SOC
    and the ΔSCF.

**Fix:**
- *Text plus re-runs on existing data:*
  - State the S1 relaxation discrepancy (0.153 against 0.076 eV) and give the adiabatic gap as a range (−0.01 to
    +0.08 eV), with the rates at both ends.
  - Use the vertical gap at the DFT S0 geometry (+0.116 eV) in the vertical-to-adiabatic step.
  - Report the ΔSCF-minimum projection for aza-BDP-Br2 (λ = 0.40 against 0.08 eV) next to the BDP and I2-BDP
    numbers. Say that the harmonic model fails along the 190–260 cm⁻¹ coordinate.
- *Cheap new calculation:* NTOs, oscillator strengths, or the overlap of the TDA roots at the DFT S0 geometry with
  the pipeline RPA S1/S2 and TDA T2. If possible, one TDA T2 geometry optimisation, or a constrained relaxed scan,
  which run_vertical_gradients.py already nearly provides.

### 2. The near-resonance tolerance is taken from proxy modes and an incorrect general argument; the candidate's own line shape gives 0.008–0.021 eV (major, text plus re-run)

**Location:** main.tex:387–388 and 410–412; tolerance_multimode.py:1–9, 54 (docstring: "since no aza-BDP-Br2 Hessian
is available"); tolerance_multimode.md.

**What is wrong:**
- *Stale proxy.* tolerance_multimode.py evaluates the aza-BDP-Br2 S1→T2 tolerance with the I2-BDP S1↔T1
  Huang–Rhys factors (λ ≈ 0.19 eV), "since no aza-BDP-Br2 Hessian is available". That Hessian and the S1/T2
  gradients now exist and give a mutual λ of only 0.014 eV. The main text (l.387–388) still quotes the proxy result,
  "0.017–0.020 eV, the same as the single-mode 0.016 eV".
- *Incorrect general argument.* l.410–412 says: "For an uphill channel the rate falls as e^{−ΔE/k_BT}, so this
  tolerance approaches k_BT ln2 = 0.018 eV whatever the line shape." That is not a general property of an uphill
  rate.
  - For a Gaussian (classical) line the logarithmic slope is (ΔE + λ)/(2λk_BT).
  - It equals 1/k_BT only when ΔE ≈ λ. For the MLJ parameters (ΔE = 0.107, λ_s = 0.10 eV) that happens by
    coincidence.
  - When λ < ΔE the slope is steeper and the tolerance tighter. With λ ≈ 0.014 eV (plus whatever λ_out is added) the classical slope is up to 3–4 times steeper; quantum hot bands soften this, and the net tolerance computed with the authors' code is 0.008–0.021 eV (table below), i.e. up to a factor of 2 tighter than 0.016–0.018 eV.
- *The resolved-comb artefact is back.* The gas-phase density for this channel is again a resolved comb, the
  artefact that round 4 raised for the single-mode model. The classical width is σ = √(2·0.0039·k_BT) ≈ 0.014 eV,
  and the rate is non-monotonic in the gap: k(0.118 eV) = 2.6e8 < k(0.156 eV) = 2.8e8 s⁻¹ (my evaluation with
  rho_multimode).

**Evidence:** tolerance_multimode.tolerance applied to the actual aza-BDP-Br2 S1/T2 Huang–Rhys factors (from
fcwd_table.aza_t2_rows):

| λ_out (eV) | Factor-2 tolerance at ΔE = +0.107 (eV) | at ΔE = +0.069 (eV) |
|---|---|---|
| 0.00 | 0.016 | 0.010 |
| 0.02 | 0.008 | 0.009 |
| 0.05 | 0.012 | 0.015 |
| 0.10 | 0.017 | 0.021 |

**Fix (text, plus re-running tolerance_multimode.py with the aza-BDP-Br2 modes):**
- Replace the proxy row.
- Replace "whatever the line shape" with the correct statement: about k_BT ln2 only when the line width is
  comparable to the gap, and tighter, down to about 0.01 eV, for the narrow, nearly parallel line computed here.
- Propagate the change to the abstract ("0.016 eV, close to k_BT ln2") and to the quadrature comparison of the
  fault-tolerant-resources section (main.tex:582–585), where 0.027 eV is then 1.3–3 times the tolerance, not 1.7.

### 3. "The class of the candidate holds (Φ_T = 0.80–1.00) and only its rate is fragile, by a factor of about 50" contradicts SI Table S15b and holds only one factor at a time (major, text)

**Location:** main.tex:416–420; si_fcwd.tex:41–52; Table S10.

**Evidence** (k_r = 2.10e8 and k_ic = 1.5e8 s⁻¹ from candidates_round1/azaBDP_Br2.json; S1→T1 channel 4e7 s⁻¹;
S1→T2 from rho_multimode with the authors' aza-BDP-Br2 Huang–Rhys factors, V2 = 46.77 cm⁻¹):

| Gap treatment (T2–S1, eV) | Gas phase k_T2 (Φ_T) | λ_out = 0.05 | λ_out = 0.10 |
|---|---|---|---|
| Pipeline vertical, +0.107 (SI Table S15b, left column) | 4.4e8 (**0.57**) | 2.2e10 (0.98) | 2.6e10 (0.99) |
| Pipeline-style vertical at DFT S0, +0.116 | 2.8e8 (**0.47**) | 1.3e10 (0.97) | 1.8e10 (0.98) |
| No T1 shift (Table S10), vertical, +0.156 | 2.8e8 (**0.47**) | 9.9e8 (0.74) | 3.1e9 (0.90) |
| No T1 shift, adiabatic, +0.118 | 2.6e8 (**0.46**) | 1.2e10 (0.97) | 1.7e10 (0.98) |
| Paper's adiabatic, +0.069 | 3.2e9 (0.90) | 1.5e11 (1.00) | 1.0e11 (1.00) |
| ΔSCF-S1 relaxation (finding 1), −0.008 | 8.1e12 (1.00) | 1.6e12 (1.00) | 6.9e11 (1.00) |

- The SI's own vertical gas-phase entry, 4.4e8 s⁻¹, already gives Φ_T ≈ 0.57. That is outside the "0.80–1.00"
  of the main text and only 0.07 above the class boundary.
- Combining two of the paper's own controls flips the class: the no-shift gap with the gas-phase multimode line
  shape gives Φ_T = 0.46–0.47. That is the outcome round 4 (4c) warned about.
- The rate spans 2.6e8 to about 8e12 s⁻¹ across these treatments, or 4.4e8 to 1.5e11 s⁻¹ even over the paper's own
  table. That is a factor of 300 or more, not "about 50".
- The SI's "Φ_T = 0.90–1.00" (si_fcwd.tex:52) is correct for the adiabatic column alone. The main text then turns
  it into a claim about the class across all treatments.
- In addition, the gas-phase entries come from the Gaussian tail of a 0.014 eV-wide line, set by a handful of
  harmonic modes below 200 cm⁻¹, which finding 1 shows to be unreliable for this dye. The gas-phase rate is
  therefore not known even to an order of magnitude.

**Fix (text-only):**
- State that the class holds in solution (λ_out ≥ 0.05 eV) for single-factor perturbations, but that it reaches or
  crosses Φ_T = 0.5 in the gas-phase multimode limit at the vertical or no-shift gap.
- Give the rate range as at least 300-fold.
- Keep the honest conclusion that aza-BDP-Br2 is the least certain prediction. Do not imply its class is robust.

### 4. "Parameter-free" multimode validation is overstated; the gradient variant used for T2 misses the measured rate by more (minor, text)

**Location:** main.tex:44–47, 322–330, 726–728; si_fcwd.tex:21–37.

**What is wrong:**
- *λ_out is a parameter.* The Conclusions (l.727) say a "parameter-free multimode model reproduces [the line
  shape] to within a factor of 2 in the rate". The outer-sphere λ_out (0–0.10 eV) is a free parameter that moves
  k(I2-BDP) by a factor of 2.8 (4.1e9 to 1.1e10).
- *Only the minima-based variant agrees.* The factor-2 agreement holds for the variant built from ΔSCF minima. The
  vertical-gradient variant, the one used for the decisive aza-BDP-Br2 T2 channel, gives 1.5e9 s⁻¹ in the gas
  phase: 4.7 times below the pipeline and 2.3 times below experiment (Table S15).
- *Inconsistent gap input.* SI S12 says that the adiabatic gap is "the consistent input of the displaced-oscillator
  line shape" (Table S15b caption), yet the S1→T1 multimode rates of BDP and I2-BDP (Table S15) use the vertical
  xTB gaps. The effect is small (−0.01 and +0.03 eV), but the convention should be the same throughout.

**Fix:**
- Drop "parameter-free" or qualify it as "parameter-free apart from λ_out".
- State that the gradient variant agrees only within a factor of about 5.
- Use the adiabatic gaps (0.75 and 0.60 eV, which exist) in Table S15, or explain the choice.

### 5. The abstract's multimode budget range omits the BDP value (minor, text)

**Location:** main.tex:44.

**What is wrong:** "0.05–0.10 eV with a multimode line shape". tolerance_multimode.json gives 0.036 eV for BDP in
the gas phase, and the main text (l.385–387) says the budget "halves in the gas phase".

**Fix:** Write "0.04–0.10 eV", or name the dye.

### 6. Archive and provenance of the multimode results (minor, reproducibility)

**Location:** paper/analysis/fcwd_multimode.json and .md; fcwd_table.py; tolerance_multimode.py.

**What is wrong:**
- fcwd_multimode.json currently holds only aza-BDP-Br2 at λ_out = 0.10 eV (`lambda_solvent_eV: 0.1`). The
  BDP/I2-BDP runs that give the 4.1e9 s⁻¹ quoted in the main text (l.326) are not archived, and the .md table does
  not state λ_out.
- The aza-BDP-Br2 entries of this JSON (k = 2.6e10 s⁻¹ at the vertical gap) are the λ_out = 0.10 values, not the
  gas-phase ones. A reader comparing them with the gas-phase discussion will be misled.
- Table S14/S15 values are regenerated in memory by fcwd_table.py. I reproduced every one to the printed digits.
- tolerance_multimode.py still says the aza-BDP-Br2 Hessian does not exist (finding 2).

**Fix:** Write one JSON per (dye, λ_out, model), including the aza-BDP-Br2 T2 rows and the ΔSCF-minimum projections,
and record λ_out in the .md.

### 7. Carry-over text items that are still wrong (minor, text)

- **l.390:** "the two SOC operators of Section 2.1 differ by 68% in V1 of I2-BDP". By the paper's own analysis
  (l.145–151) the difference is between *states* (ECP/RPA against all-electron/TDA); the operators agree within 7%.
  Write "the two SOC treatments".
- **l.374–375:** "up to a factor of 2.7" needs "except aza-BDP-Br2 (×77)" (round-4 #8).
- **l.358–361 and 424–425:** the I2-BDP-0 alternatives that come "close to experiment" (Φ_T = 0.67 and 0.73) use
  V2 = 8.0 cm⁻¹, which the all-electron benchmark puts at 3.2–3.5 cm⁻¹ (Table S6). With the benchmark coupling both
  give about 2e8 s⁻¹ (Φ_T ≈ 0.3) (round-4 #9). This also weakens the argument (l.421–425) that the S1-shift model
  would rescue I2-BDP-0.
- **l.314–315 and 709:** the functional dependence is still presented as a gap effect. PBE0 and BHandHLYP lower V1
  to 10.2 and 2.8 cm⁻¹ (validation/runs_xc). Add V1 to Table S12 (round-4 #13).
- **l.192:** the δ notation (round-4 #14).
- **si.tex:138, 142:** the (14,14) rows (round-4 #12).
- **si.tex:174, 200:** the decision cost still assumes a 0.5 eV S1–S2 gap. For aza-BDP-Br2 it is 0.10 eV (RPA) or
  0.13 eV (TDA at the DFT geometry), where the paper's own sensitivity row gives C_S1 = 1.3e8 Toffolis. This changes
  the aza-BDP-Br2 cost per call by about 50% (round-4 #4e).
- **Table S6 (si.tex:259–272):** the T5 column of the "eff." row for I2-BDP reads 1119 cm⁻¹ (and SOMF 848).
  validation/runs/I2-BDP.json confirms V5 = 1119.4 cm⁻¹. A ππ*–T5 coupling of 1100 cm⁻¹ is extraordinary: T5 must
  be an iodine-centred n/σ* state. Say so in the caption, since it is the largest coupling in the paper and the
  reason T5 is excluded from the sum.

### 8. Still missing for the central claims (major, carried over; needs calculation)

The two gaps on which the narrative turns, T2–S1 of I2-BDP-0 and of aza-BDP-Br2, still have no correlated
single-reference reference. The bromine couplings of the top candidate still have no all-electron benchmark, so its
V2 = 46.8 cm⁻¹, the largest S1–T(2,3) coupling in the set at ζ(Br)/ζ(I) ≈ 0.5, is unchecked. Both are listed only as
limitations (l.716–717). Given findings 1–3, a single SOMF calculation on aza-BDP-Br2 at the pipeline geometry
(run_soc_benchmark.py already exists) and one ADC(2)/SCS-ADC(2) or STEOM-CCSD calculation of S1/T1/T2 for that dye
would decide whether the "top candidate" survives. I regard the SOMF check as required and the ADC(2) as strongly
recommended.

---

## Numbers re-verified this round (no action needed)

- Table S11: every entry reproduces from tda_full_check.json (ΔT_n − ΔT_1).
- Table S9 against decision_cost.json: all D1 and D2 values.
- Table S14/S15 (BDP, I2-BDP; minima and gradients; λ_out = 0, 0.05, 0.10 eV): reproduced with
  fcwd_multimode.analyse. I2-BDP k_multi = 4.10e9, 7.96e9 and 1.14e10 s⁻¹; λ_S1T1 = 0.129 and 0.186 eV;
  ΣS = 1.15.
- aza-BDP-Br2 S1→T2 table (si_fcwd.tex:44–50): 4.4e8/3.1e9, 2.2e10/1.5e11, 2.6e10/9.9e10 s⁻¹; λ = 0.153, 0.086 and
  0.191 eV; mutual λ = 0.014 eV; MLJ 1.2e10 s⁻¹.
- Multimode tolerances (l.384–388) match tolerance_multimode.json.
- Adiabatic aza-BDP-Br2: S1 vertical-to-adiabatic 0.08 eV, relaxation energies 0.076/0.094 eV, S1–T1 0.68 to
  0.70 eV.
- Table S10 aza-BDP-Br2: no-shift k = 1.44e9 s⁻¹ (×8.7), Φ_T = 0.80.
- SOC text: 22.2/13.2 → ×2.8; 8.0/3.5 → ×5.2; summed V² 1135/1108 (2%) and 844/1007 (16–19%).
- Quadrature 0.027 eV = 1.7 × 0.016 eV (as written; see finding 2 for the tolerance itself).
