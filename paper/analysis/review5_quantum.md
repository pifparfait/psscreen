# Referee report (quantum algorithms and resource estimation), round 5

Manuscript: "How Accurate Must Excited States Be to Screen Triplet Photosensitizers?" (main.tex 2026-11-25 19:23,
main.pdf and si.pdf rebuilt 19:23). Line numbers refer to main.tex and si.tex as they stand now. si.tex is
byte-identical to the round-4 version. si_tables.tex changed only in the Table S7 caption. No script in paper/quantum/
and no file in paper/quantum/results/ changed after 2026-11-22, so every new number in the text has to follow from the
round-4 outputs. I checked each one against decision_cost.json/md, alpha_so.json and *.cas_p.json. My checks ran in
/tmp/claude-0/review5q/.

## Summary

The revision adds a classical-baseline and active-space caveat to the abstract, the MLAE coherent depth, the cost of
measuring gamma^2, a regenerated Table S9, and narrower claims for the D2 test and the MLAE advantage. Most of the new
numbers re-derive correctly: 1/200, 0.4-50 days, 4-140x, 1e3 calls, 3e10 Toffolis, the CASCI timings and Table S9.
Three of the new statements are still wrong or incomplete:
- the gamma^2 paragraph understates its own consequences (its cost share, the fact that D2 needs gamma^2 too, and
  the error it adds to the 20%);
- the "within a factor of 2" certificate is still claimed in the body and in the Table 3 footnote;
- the "40-fold" spin-orbit variability is not a like-for-like comparison.

None of round 4's new-calculation requests (findings 1b, 2b and 3b) was done; the authors replaced them with caveats.

## Recommendation

**Minor revision.** The authors chose to reframe rather than recompute: the abstract now says the costs are for a
classically solvable space and are not a quantum-advantage claim. I accept that position if it is carried through the
body and the Conclusions, and if the errors below are fixed. Every remaining item is text-only or a calculation of a
few lines. If the editor wants the paper to make any statement about cost at a classically hard, SOC-aware active
space, round-4 findings 1b and 2b still apply, and the recommendation becomes major.

## Status of round-4 findings

| # | Finding (round 4) | Status | Comment |
|---|---|---|---|
| 1 | Costed problem is classically trivial; no classical baseline | Partially resolved | Abstract l.57-59 now says this. Sec. "Cost of a screening decision" (l.587-632), Discussion and Conclusions (l.730-733) do not. No calculation at a classically hard size (1b). See N4. |
| 2 | Three inconsistent models; alpha_SO varies 40x | Partially resolved | Abstract caveat added, but the "40-fold" mixes orbitals (N3). No cost with an SOC-aware space. The per-dye lambda (12.0 and 11.5 Ha) is stated but not used (l.580). |
| 3 | Coherent depth and runtime understated | Partially resolved | MLAE deepest circuit added (l.630-632), arithmetic verified. Still missing: the canonical depth in years (4.5e14 Toffolis is 14 years), D1 wall-clock range (up to 28 years), campaign time, parallelism over MLAE depths, a rigorous low-depth estimator, physical qubits. See N5. |
| 4 | "Within a factor of 2" is not what D2 delivers | Partially resolved | Abstract fixed (l.55-56). Body l.627-628 and Table 3 footnote a (l.656-657) unchanged. See N2. |
| 5 | gamma^2 must be measured; polynomial bias not in the 20% | Partially resolved | Measurement cost added (l.236-238), with errors (N1). The 10% polynomial bias is not addressed. The CAS(4,4) exact gamma^2 = 0.94 is stated (l.550). |
| 6 | si.pdf stale; S9 caption vs Table 3 | Resolved | si.pdf matches decision_cost.json for every row. The "by MLAE" caption is consistent with the D2 values shown (aza-BDP-Ph4 1.6e11, BDP-0 1.8e11). One residual point in N7. |
| 7 | Definition of p (eq 3 vs text vs Fig. 2 axis) | Not addressed | Eq 3 (l.220) has p = k/K', l.592 and Table 3 have p = gamma^2 k/K'. The axis of figures/decision_cost.pdf still reads "p = k_ISC/K'" (decision_cost.py:489), while the Fig. 2 caption says gamma^2 k/K'. |
| 8 | Bare-Hamiltonian comparison omits the core | Resolved | l.599 now reads "1200 times smaller ... 27 and 800 times larger". Verified: 1250, 27, 815. |
| 9 | Walk steps per call understated | Resolved | l.569, "2-6e4 walk steps". |
| 10 | alpha/||chi|| = 11 uses a TDDFT chi | Not addressed | l.600 and l.703-706 unchanged. The CASCI chi over T1-T6 is 114 cm^-1 (I2-BDP.cas_p.json), which gives alpha/||chi|| = 105, not 11. |
| 11 | S1-preparation assumptions (gap, gamma^2, leakage) | Not addressed | si.tex l.172-179 unchanged. The CASCI S2-S1 gap of aza-BDP-Br2 is 0.34 eV, not 0.5. |
| 12 | Qubit and walk bookkeeping | Partially resolved | Text gives "19-20 phase bits" and 220, and states the uncounted registers (l.572-574). logical_qubits.json still has degree 51264 and 19 bits. decision_cost.md still says "fitted". The rank-33 row is still missing from thc_pennylane_resources.md. |
| 13 | C_SO counts spatial orbitals | Not addressed | Negligible for the totals. |
| 14 | MLAE constant and its range | Not addressed | No Cramer-Rao justification. The use of the constant at relative errors far above 1 for light-atom D2 is not flagged. |
| 15 | D2 sampling vs MLAE not like for like | Not addressed | SPRT is costed by its expected sample number, MLAE by a fixed quantile. The stage schedule uses L from the true theta. |
| 16 | Missing error terms in the near-resonance budget | Partially resolved | Quadrature updated with the new TDA value; 0.027 eV verified. Block-encoding precision, dependence of d on lambda_s, and the 95% factor for QPE are still missing. |
| 17 | Channel route: favourable assumptions | Not addressed | |
| 18 | Abstract optimistic without saying so | Not addressed | "Calls of about 2e8" (actually 0.9-2.2e8). The abstract does not say the totals are MLAE or best-estimator values. |
| 19 | Missing prior work | Not addressed | None of the suggested references is cited. Fomichev2024XAS is already in related_work.bib but uncited. |

## New findings

### N1. The gamma^2 measurement paragraph is wrong in three ways (moderate)

**Location:** main.tex:236-238; decision_cost.py:264; si.tex:172-179 (not updated).

**Problem and evidence.**
- **(a) The 0.4% bound fails for one Table 3 column.** "At most 0.4% of any D1 cost below" holds for the MLAE and
  canonical columns: 3e10/7.2e12 = 0.42% for aza-BDP-Br2. It fails for the D1 channel column of Table 3, which needs
  the same S1 normalization (its success probability is gamma^2 V_n^2/alpha^2, l.637). There gamma^2 is 3e10/4.1e11
  = 7.3% for aza-BDP-Br2, 2.2% for aza-BDP-BrPh2 and 1.8% for Cl2-BDP.
- **(b) D2 also needs gamma^2.** The text says gamma^2 is needed "since a D1 estimate needs gamma^2". But the D2 test
  compares gamma^2 p with gamma^2 p_th (decision_cost.py:264: `p, p_th = GAMMA2 * x["p"], GAMMA2 * x["p_th"]`), so
  gamma^2 is needed to place the zone edges. This cost is not included. 3e10 Toffolis is 23% of the D2 cost of
  aza-BDP-Ph4 (1.3e11), 19% for BDP-0, 10% for BDP and 8% for Cl2-BDP. These are the light-atom negatives, whose
  cost range ("1e11-1e12") the abstract quotes.
- **(c) The extra error is not in the D1 budget.** A 10% error on gamma^2 enters k = K'(gamma^2 p)/gamma^2 directly.
  - In quadrature with the 20% on gamma^2 p, the error is 22%.
  - With a union bound (joint confidence 0.90) and the deterministic polynomial bias of up to r/2 = 10% that round-4
    finding 5 pointed out, the worst case is 1.2 x 1.1/0.9 = 1.47. That is ln 1.47 = 0.39, 56% of the ln 2 budget,
    not the "about a quarter" claimed at l.226.
- **(d) Minor points.**
  - At p = 0.5, MLAE does not need the small-angle constant: dtheta = c/(2N) and dp/p = 2 dtheta give N = 450. Direct
    sampling to 10% at 95% needs only z^2(1-p)/(p r^2) = 384 shots of C_S1. "About 1e3 calls" is conservative, which
    is fine.
  - The SI cost-model section (S4) still describes gamma^2 only as an assumption.

**Fix.**
- Measure gamma^2 to 3% instead. That is about 1e4 shots, about 1e11 Toffolis, at most 1.4% of any D1 MLAE cost. Its
  contribution is then negligible in quadrature: sqrt(0.2^2 + 0.03^2) = 0.202.
- Either tighten the polynomial precision to eps = r sqrt(p)/20, which adds only about ln 5 to the degree slope, or
  state D1 as "20% statistical plus at most X% bias at confidence 1 - 2 delta".
- Add the gamma^2 cost to D2, or say that D2 needs only about 30% on gamma^2 (a few hundred shots).
- Restrict "0.4%" to the filter route.
- Mirror the paragraph in SI S4.

Text plus a small calculation, no rerun of the model needed.

### N2. The inside-zone "certificate" is still claimed in the body and in Table 3 (moderate, text-only)

**Location:** main.tex:627-628 ("the test ends ... with the answer that the rate is within a factor of 2 of the
threshold") and 656-657 (Table 3 footnote a, the same wording). These contradict l.228-229 ("may be returned in
either class") and the corrected abstract.

**Evidence.** Recomputed from decision_cost.json with the code's own rule, dth = max(theta - theta0, theta1 - theta)/2.

| Dye | dth | Near-edge exclusion needs dth < | Guaranteed? |
|---|---|---|---|
| Br2-BDP | 1.01e-4 | 5.9e-5 | No |
| aza-BDP-IPh2 | 1.82e-5 | 1.50e-5 | No |

So the stopping rule certifies one side only, for example "k < 2 k_th".

**Fix.** Replace both passages with "the test returns a class that is correct up to the factor-2 indifference zone".
Alternatively, cost the two-sided rule dth < min(.)/2 for these two dyes (a small calculation).

### N3. "Varies 40-fold" compares different orbitals and understates the range (minor, text-only)

**Location:** main.tex:58 (abstract); compare l.602-603.

**Evidence.**
- 40 = 12006 cm^-1 (canonical window, B3LYP orbitals, alpha_so.json) / 299 cm^-1 (AVAS, Hartree-Fock orbitals,
  Table 2). With the same Hartree-Fock orbitals (Table 2), the canonical/AVAS ratio is:

  | Dye | Ratio |
  |---|---|
  | I2-BDP | 28 (l.602-603 itself quotes 8300 against 300) |
  | I3-BDP | 38 |
  | I-BDP | 70 |
  | I-BDP-0 | 28 |
  | I2-BDP-0 | 1.5 |

- Changing only the orbitals in the same canonical window, from Hartree-Fock to B3LYP, already changes alpha_SO by a
  factor of 1.4 for I2-BDP (8315 to 12006) and 6.4 for I2-BDP-0 (418 to 2656). The cost of I2-BDP-0, the most
  expensive dye, therefore moves by 6x through the orbital choice alone.

**Fix.** Write "by a factor of up to about 80 (28 for I2-BDP)", and mention the HF/B3LYP sensitivity in the
paragraph at l.602-607.

### N4. The classical-baseline caveat appears only in the abstract, and "how the cost scales" is not supported (minor, text-only)

**Location:** main.tex:57-59, 587-607, 730-733.

**Evidence.**
- The CASCI timings (407, 983 and 2341 s, *.cas_p.json) support "minutes".
- The phrase does not appear in Section "Cost of a screening decision", the Discussion or the Conclusions. The
  Conclusions still state "On a fault-tolerant computer a screening decision costs 1e11-7e15 Toffolis" without
  qualification.
- "They measure how the decision cost scales" is not supported: the paper has one system size and no N-dependence.
  What it does establish is how the cost depends on alpha_SO/V_eff and on the distance to the threshold (Fig. 2,
  l.701-706).
- Exact CASCI and QPE on the same bare Hamiltonian give the same rate, 1200 times off for I2-BDP (l.596-599). The
  costed machine therefore corresponds to no Hamiltonian that either solver has in hand. This is stronger than "not a
  quantum advantage", and the reader should be told.

**Fix.**
- Replace "how the decision cost scales" with "how the decision cost depends on the spin-orbit normalization and the
  distance to the threshold".
- Repeat the caveat in one sentence at l.595 and in the Conclusions.

### N5. MLAE coherent depth verified, but the runtime picture is still incomplete (minor, text-only)

**Location:** main.tex:628-632.

**Verified.** The ae_demo.py schedule is m = 0, 1, 2, ..., 2^K with 100 shots per depth. Total calls are about
100 x 2^(K+2) and the deepest circuit is 2^(K+1) + 1 calls, so the ratio is 1/201 at K = 8. D1 MLAE divided by 200
gives 3.6e10 Toffolis (aza-BDP-Br2, 0.42 days) to 4.4e12 Toffolis (I2-BDP-0, 50 days). These match the text.

**Missing.**
- **(a) Logical error rate.** A 4e12-Toffoli circuit needs a logical failure rate of about 1e-13 per Toffoli. The
  text gives 1e-15 for canonical estimation but implies that MLAE is benign.
- **(b) Wall-clock time.** "About 6 years" is quoted for I2-BDP only. The D1 MLAE range is 0.2 years (aza-BDP-Br2)
  to 28 years (I2-BDP-0), and 14.5 years for aza-BDP-IPh2. The deepest canonical circuit of I2-BDP-0 is 14 years of
  coherent time. The D2 campaign over the 15 dyes totals 1.2e15 Toffolis, 38 years serially.
- **(c) Parallelism.** MLAE shots at different depths are independent, so the 200-fold depth reduction also allows
  parallel execution. Say so, since it is the main operational advantage.

**Fix.** Add one sentence on each point.

### N6. The abstract and the body give different MLAE-advantage ranges, and the light-atom range is clipped (minor, text-only)

**Location:** main.tex:53-57, 619-622.

**Evidence.**
- The D2 sampling/MLAE ratios are:

  | Dye | Ratio |
  |---|---|
  | aza-BDP-Br2 | 3.6 |
  | aza-BDP-BrPh2 | 4.5 |
  | Br2-BDP | 36 |
  | I1-BDP | 61 |
  | I2-BDP, I3-BDP | 69, 71 |
  | I2-BDP-0 | 102 |
  | aza-BDP-IPh2 | 142 |

  So the abstract's "4-140" is right for the ISC-active Br/I dyes. But the abstract does not say that this is the
  class decision; for D1 the ratio is 3e2-2.5e4. The body (l.621) says "4-270" for "heavy-atom dyes", which includes
  I-BDP (274).
- The abstract gives "1e11-1e12 for light-atom negatives", while the body gives 1e11-2e12; the BODIPY core is 2.2e12.
- "Calls of about 2e8" should be 0.9-2.2e8 (round-4 finding 18).

**Fix.** Write "for the class decision" and use one dye set in both places. Write 1e11-2e12 and 0.9-2.2e8.

### N7. The SI does not reflect the revision (minor, text-only)

**Location:** si.tex (unchanged since round 4); Table S9.

**Problem.**
- Section S4, to which the main text sends readers for "details and sensitivities" (l.249), has neither the gamma^2
  measurement nor the MLAE depth schedule.
- Table S9 lists D1 costs for the dyes with k < 1e6 s^-1, up to 8.2e18 (MLAE) and 4.7e19 (canonical) for I-BDP. The
  main text says D1 is evaluated only for k >= 1e6 and quotes totals up to 7e15.

**Fix.** Add the two paragraphs to Section S4, and put dashes or a footnote ("not a relevant decision") in the
Table S9 D1 cells of the slow dyes.

### N8. Cross-check outside my remit: the abstract's multimode tolerance range (minor, text-only)

**Location:** main.tex:44 against l.384-387.

**Problem.** "0.05-0.10 eV with a multimode line shape" covers the I2-BDP values (0.046-0.098 eV) but not BDP
(0.036 eV in the gas phase, 0.068-0.069 eV with an outer-sphere term). The body says the 0.08 eV "halves in the gas
phase".

**Fix.** Write "0.04-0.10 eV", or "for I2-BDP".

## Checks that passed

- **CASCI timings and p_CAS ratios.** 407, 983 and 2341 s. The p ratios are 1250 (I2-BDP), 27 (aza-BDP-Br2) and
  815 (core).
- **MLAE coherent depth.** 1/200 of the calls; 3.6e10 to 4.4e12 Toffolis; 0.42 to 50 days.
- **Cost of gamma^2 by MLAE.** 636 calls x C_S1 = 1.7e10 Toffolis, about 3e10 as stated. 0.42% of the cheapest D1
  MLAE cost.
- **D2 ratios.** Sampling/MLAE is 3.6-142 (ISC-active Br/I dyes) and 3.6-274 (all Br/I dyes). Canonical/sampling is
  1.24-9.4 for the six dyes named. The cost of an unknown margin is 2.2-5.3.
- **Table 3 against decision_cost.json.** Every entry matches. Table S9 in si.pdf matches the json, row by row.
- **Section 5 ranges.** D1 MLAE 7.2e12-8.7e14, canonical 3.9e13-7.3e15, 4.8e4-4.0e6 calls, channel route 4.1e11-7.5e13
  (filter/channel ratio 11.6-28.8), alpha/V_eff 150-1.27e4, Zhou et al. ratios 190-1.9e5 (D1) and 3.4-58 (D2).
- **Near-resonance quadrature.** sqrt(0.003^2 + 0.027^2 + 0.003^2) = 0.027 eV, 1.7 times the tolerance.
- **Qubit totals.** logical_qubits.json gives 201 (MLAE) and 220 (canonical, 19 bits). With 20 bits the canonical
  total is 221, which is still "220".
