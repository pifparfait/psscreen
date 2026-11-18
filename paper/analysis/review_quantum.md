# Referee report (quantum algorithms and resource estimation): "How Accurate Must Excited States Be to Screen Triplet Photosensitizers? ..." (JCTC draft)

Scope of this report: the gap-resolved ISC observable (Methods sec:filter-method, Results sec:filter,
Table `tab:filter`), the decision-level cost model (sec:decision-method, eq `eq:p`, sec:resources,
sec:decision, Table `tab:decision`, Figure `fig:decision`, SI Figure S3), the THC resource numbers, and the
comparison with Zhou et al. (PRA 2026 / arXiv:2512.15889) and Eklund et al. (arXiv:2603.19007). Line numbers
refer to the current `paper/main.tex` (726 lines). I checked every number against `paper/quantum/*.py`,
`paper/quantum/results/*` and `psscreen/isc_filtered.py`, and re-ran one part of `degree_check.py` (see M7).
Points from the previous report (`analysis/manuscript_review.md`) that the text now fixes are not repeated.
Points that remain open are marked **[open from previous report]**.

## Summary and recommendation

The quantum part makes a useful and, to my knowledge, new point for this application. It writes the
golden-rule ISC rate as a QSVT spectral filter over all triplets and checks it in simulation. It then costs the
screening decision rather than the shot, and shows that the success probability is set by the spin-orbit
block-encoding norm and is very small (6×10⁻⁹–4×10⁻⁵). So the number of calls dominates, and amplitude
estimation (AE) is needed. The arithmetic is reproducible. Every entry of Table `tab:decision` and every
range in sec:decision matches `decision_cost.md`, except the few listed under Minor points.

The cost model does have several logical problems:

- The S1 state preparation is not what the text says it is.
- The "cost lever" in the Discussion uses the wrong quantity.
- The confidence and family-wise statements rest on an empirical MLAE constant outside the range where it
  was measured.
- The S1 filter ignores near-degenerate triplets.
- The QPE precision uses up most of the accuracy budget.
- An obvious alternative (channel-resolved energies and couplings with the golden rule evaluated classically)
  is not compared.

None of this overturns the qualitative message, "calls, not cost per call, dominate". Several items move the
headline numbers by factors of 2–10. Two of them change what the paper should recommend.

**Recommendation: major revision** (quantum part).

---

## Major points

### Q1. The S1 preparation: one round of amplitude amplification at γ² = 0.5 does nothing, and the 1/γ scaling claimed in the text is not in the model

**Where:** l. 229–231, l. 241–243; `decision_cost.py` `amp_rounds`, `c_s1`; Limitations l. 676.

**Problem.**
- **The amplification gives no gain.** With γ² = 0.5, θ = asin γ = π/4. After m = 1 round the success
  amplitude is sin(3θ) = sin(3π/4) = 0.707, which is the same as with no amplification. The optimal real m is
  π/(4θ) − ½ = 0.5, and `amp_rounds` rounds it up to 1.
- **The model charges for it anyway.** It pays 3× the preparation cost (C_S1 = 6.6×10⁷, 33–55% of one call)
  for no gain. It then treats the preparation as deterministic: p in eq `eq:p` carries no factor for the S1
  filter flag.
- **The true success probability is about half the stated one.** An eigenstate filter applied to a
  sum-of-Slaters state with overlap γ² succeeds with probability γ² (here 0.5 even after the "amplification").
  The good-state probability that AE sees is about γ²·p, not p.
- **The 1/γ statement is not what the code computes.** "They scale … through the amplitude amplification of
  the S1 preparation, as 1/γ" (l. 242–243) does not describe the implementation. In the code γ enters only
  through a step function in `amp_rounds` and as an additive cost.

**Fix.** Choose one consistent design and cost it.
- (a) Fold the S1 filter flag into the AE good state: p → γ² p, the number of calls grows by 1/γ, and C_S1 is
  paid once per call. This is the design that has the 1/γ scaling, and probably the cheapest.
- (b) Alternatively, use fixed-point amplitude amplification (Yoder–Low–Chuang), with cost
  O(log(1/δ)/γ) applications, or exact AA with a tuned final phase (this needs γ to be known).

With (a) the per-call cost falls from about 2.0×10⁸ to 1.6×10⁸ and the number of calls rises by √2, so the
totals change by about +10%. The logic, though, is currently wrong. The sentence at l. 242–243 must also
describe what is actually computed, and the γ sensitivity (γ² = 0.1, 0.5, 0.9) should be tabulated.

### Q2. The S1 eigenstate filter must resolve S1 from nearby triplets, not only from S2

**Where:** l. 229–231 ("an S1 eigenstate filter for a 0.5 eV S1–S2 gap"); `decision_cost.py` `DELTA21_HA`.

**Problem.**
- **The filter cannot tell singlets from triplets.** It is a polynomial in the spin-free H. It passes every
  eigenstate within its width of E(S1), including M_S = 0 triplet components.
- **The trial state will contain such components.** A sum of 10⁴ Slater determinants that is not
  spin-adapted generically has M_S = 0 triplet admixture.
- **The nearby triplets are exactly the ones that matter.** The triplets closest to S1 drive ISC: 0.11 eV for
  azaBDP_Br2, 0.06–0.25 eV in the CAS simulations of `tab:filter`, and 0.17 eV for aza-BDP-Ph₄ at NEVPT2.
- **Leakage biases the rate directly.** A T_n component that leaks into the prepared "S1" passes through H_SO
  and the ρ filter, and biases k_ISC.

With a filter gap of 0.1 eV instead of 0.5 eV, the S1-filter degree (and C_S1) grows 5×. C_S1 then becomes
about 3×10⁸, the largest term in a call.

**Fix.** Do one of the following and cost it: (i) prepare a spin-adapted (CSF) trial state; (ii) add a spin
projection, for example a block encoding of Ŝ² with its own filter; or (iii) set the filter gap to
min(E_S2 − E_S1, |E_Tn − E_S1|) for every T_n with non-zero overlap. In each case, state the residual triplet
leakage as an error in k_ISC.

### Q3. The "cost lever" in the Discussion uses ‖χ‖, but ‖χ‖ does not set p

**Where:** l. 656–660 ("p ∝ (‖χ‖/α_SO)², the decision cost scales as α_SO/‖χ‖ … would reduce the cost by the
same factor"); l. 558–561 and the ‖χ‖ column of `tab:decision`.

**Problem.**
- **What sets p.** By eq `eq:p`, p = k_ISC/K′ = (V_eff/α_SO)², with V_eff² = Σ_n V_n² ρ(ΔE_n)/ρ_max, the
  Franck–Condon-weighted coupling. ‖χ‖ enters only as an upper bound (p ≤ ‖χ‖²/α_SO²).
- **‖χ‖ is not what makes p small.** From `decision_cost.json`, V_eff = √p·α_SO is 4.2 cm⁻¹ for I2-BDP, against
  ‖χ‖ = 1120 cm⁻¹. So 1/√p = 2850 splits into α_SO/‖χ‖ = 11 times ‖χ‖/V_eff = 270. For I2-BDP-0 the split is
  65 × 195, and for azaBDP_IPh2 156 × 44. The dominant suppression is that most of H_SO|S1⟩ goes to triplets far
  from resonance, which the filter rejects.
- **Consequences for the Discussion.** The lever is α_SO/V_eff, not α_SO/‖χ‖, and the factor "5–160" is not
  the possible gain. An encoding restricted to the near-resonant couplings could gain up to α_SO/V_eff
  (150–1.3×10⁴). A "resonance-weighted preparation of χ" cannot remove the ‖χ‖/V_eff part without
  pre-filtering χ, which is the filter itself.
- **Mixed levels of theory [open from previous report, M5(c)].** ‖χ‖ comes from TDDFT over a few computed
  triplets. α_SO comes from the CAS(12,11) window. The ratio α_SO/‖χ‖ is therefore still a mixed-level
  quantity.

**Fix.** Rewrite the paragraph at l. 656–660 in terms of α_SO/V_eff (equivalently 1/√p), and give the split
α_SO/‖χ‖ × ‖χ‖/V_eff per dye. Replace the ‖χ‖ column of `tab:decision` by V_eff, or add V_eff. Compute ‖χ‖ in
the same CAS space as α_SO, or drop it.

### Q4. Missing baseline: channel-resolved energies and couplings, with the golden rule evaluated classically

**Where:** sec:filter-method, sec:decision; Discussion l. 640–654.

**Problem.** In these dyes 1–3 channels carry the rate (Table `tab:isc`, the budget analysis). An obvious
alternative estimates E(T_n) by QPE (1.1×10⁸ Toffolis each at 2 mHa, the paper's own number) and V_n² by AE on
an eigenstate-filtered H_SO|S1⟩. The golden rule is then evaluated classically.
- **More probability per call.** For the T_n channel the success probability is V_n²/α_SO², not
  V_eff²/α_SO². For I2-BDP that is (13.2/12006)² = 1.2×10⁻⁶, against p = 1.2×10⁻⁷, so about 3× fewer calls.
- **A cheaper filter.** The eigenstate filter needs width |E_Tn − E_Tn±1| (tenths of an eV), not the MLJ width
  of 0.07 eV, so its degree is several times lower.
- **The line shape stays classical.** The route keeps the MLJ parameters (λ_s, ħω, S) classical and
  re-evaluable. That matters because the paper's own budget shows a 0.03–0.13 eV dependence on λ_s
  (Fig. `fig:tolcurve`b).
- **It gives what the budget asks for.** The individual gaps are the quantity the accuracy budget is about.

As written, the paper's argument for the filter ("absolute, mechanism-resolved rate", l. 598) is made only
against the gap-free proxy, not against this route. The filter wins only when many triplets contribute
comparably.

**Fix.** Add a short cost comparison of the channel-resolved route for I2-BDP, azaBDP_Br2 and Br2-BDP, using
the same model as `decision_cost.py`. State when the filter is preferable: a dense triplet manifold, or
unknown channels.

### Q5. Confidence levels: the empirical MLAE constant is used outside the range where it was measured, and the family-wise factor is unsupported for MLAE

**Where:** l. 234–239, 570–581, 592–593, 676–677; SI Fig. S3; `ae_demo.py`, `decision_cost.py`
`mlae_calls`.

**Problem.**
1. **The constant comes from a simulation with a floor.** c = 45 is the median of c_eff over depths K = 7–13
   at a single p = 5.6×10⁻⁶ (`ae_demo.json`). At the true amplitude (`q95_true`) the 95th-percentile error does
   not fall below 6–9% for K ≥ 10, and it rises again from 7.2% at K = 9 to 9.4% at K = 13. Only the ensemble
   reaches 3%. The text says "stops improving at 3% … used only for errors of 5% and above" (l. 580–581). For
   the actual amplitude the floor is above 5%. Br2-BDP D2 needs about 3% angle precision (a 6% margin in p),
   which lies inside the failing range.
2. **The floor is probably an artefact of this estimator.** It is not an intrinsic feature of MLAE. The
   likely cause is the hierarchical grid maximum-likelihood search, whose window w = 30/√Fisher can lock onto
   a wrong branch. The cause is still not stated **[open from previous report, minor 3]**.
3. **The family-wise factor is not supported.** The factor ln(40N)/ln 40 (l. 592–593) assumes a
   log(1/δ)-type tail. For MLAE that tail has not been shown. The simulation uses 400 trials (16 amplitudes × 25
   for the ensemble), which cannot resolve the 99.95th percentile that N = 100 at family-wise 0.05 needs, and
   the observed q95/median ratios of 3–15 at large K point to heavy tails.
4. **D2 uses knowledge of the answer.** For D2 the number of calls is set from the true margin
   |θ_th − θ|, which is unknown before the measurement. SPRT is also run against the true p. A real test is
   sequential or adaptive, with a log-factor overhead.
5. **The canonical-QAE constant is not used.** The text says the constants "45 for MLAE and 109 for canonical
   estimation are used in the cost model" (l. 579–580). The code uses only 45. Canonical QAE is costed from
   the Brassard bound (Thm 12), which gives 2.7× more calls than the empirical 109. The canonical totals are
   never shown in the paper. They are 6–8× the MLAE totals: D1 I2-BDP 8.2×10¹⁴, D2 Br2-BDP 1.7×10¹⁵, D2
   azaBDP_IPh2 3.5×10¹⁵ (`decision_cost.json`). With the rigorous estimator the upper end of the headline
   becomes 10¹⁵–4×10¹⁵.

**Fix.**
- Use an AE variant with proven confidence intervals and no QPE, for example IQAE (Grinko2021, already cited)
  or the power-law/Rall–Fuller schemes. Cost with its rigorous constant, or at least report the canonical totals
  next to MLAE in `tab:decision` and in the abstract range.
- Fix or explain the MLAE floor.
- Derive the family-wise factor for the estimator actually used.
- State that the D2 costs assume an oracle margin, and give the overhead of an adaptive version.
- Correct l. 579–580.

### Q6. The QPE that centres the filter uses up most of the accuracy budget

**Where:** l. 231–232; `EPS_E_HA = 2e-3`.

**Problem.**
- **The QPE alone uses most of the budget.** The budget allows about 2 mHa per state, 0.06–0.08 eV on the
  gap (l. 376). The model gives the QPE of E(S1) the same 2 mHa (0.054 eV) and calls the result a systematic gap
  error (l. 232). That one error then uses 65–90% of the factor-2 tolerance, before any Hamiltonian, THC or
  active-space error.
- **For the top candidate it is over budget.** For azaBDP_Br2 (tolerance 0.016 eV, l. 392) it is 3.4× the
  budget on its own.

**Fix.** Allocate a stated fraction of the budget to the QPE centring, for example ε_QPE = 0.3–0.5 mHa. This
costs 4–7× the QPE (4–8×10⁸ Toffolis once per dye), which is negligible against the totals. Write the full error
budget as statistical (r) + QPE centring + THC/precision + Hamiltonian, so that the parts sum to the factor-2
tolerance.

### Q7. Unvalidated THC and filter-degree inputs, and one λ for all dyes

**Where:** l. 225–228, 257–260, 541–553; `thc_pennylane_resources.py`; `degree_check.py`;
`decision_cost.py` (`LAMBDA = 15.9`, `SHIFT` from the BODIPY core).

**Problem.**
1. **No THC factorization was run.** `qre.THCHamiltonian(num_orbitals=11, tensor_rank=22, one_norm=15.9)`
   takes λ from double factorization ("the THC one-norm is of the same size", docstring). Neither the THC λ
   nor the error at rank 22 or 44 has been computed. Rank 2N is low for small active spaces. The error on the
   S1–T_n gaps is still missing from the budget **[open from previous report, M8; TODO at l. 552]**.
2. **One λ′ for every dye.** λ′ = 19.4 Ha is the BODIPY-core value and is used for all 15 dyes. The canonical
   windows of the iodinated dyes contain iodine lone pairs and σ* orbitals (Table `tab:window`), so their λ
   (and λ′) differ. The filter degree is linear in λ′.
3. **The degree fit is extrapolated.** It uses three precisions (10⁻², 2×10⁻³, 10⁻⁴) and λ/Δ ≤ 1600. It is
   applied at ε = 4×10⁻⁶–2×10⁻⁵ and λ′/Δ = 7340. I evaluated `degree_check.min_degree` directly:
   - λ/Δ = 1600, ε = 10⁻⁵: degree 11 028 (fit gives 10 035, so the fit is 10% low);
   - λ/Δ = 1600, ε = 3×10⁻⁶: degree 11 868 (fit 11 103, 7% low);
   - λ/Δ = 400, ε = 10⁻⁵: degree 2756 (fit 2509, 10% low).

   The effect is small, but the fit underestimates in exactly the regime that is used. The slope is also
   neither linear nor √-linear in ln(1/ε) (slope/√ln = 1.23, 1.22, 1.68). Note too that the `degree_check.py`
   docstring still describes the old preprint-based model.

**Fix.** Run an actual THC factorization (with BLISS if Zhou et al. used it) at rank 22 and 44 for the BODIPY
core and I2-BDP. Report λ_THC and the gap errors, and use the per-dye λ′. Compute the degree directly at the
ε and λ′/Δ that are used (about 5×10⁴, which is cheap) instead of extrapolating.

### Q8. Wall-clock time, coherent depth and logical qubits are not reported, and they qualify "required in every case"

**Where:** abstract l. 53–57; sec:decision; Conclusions l. 691–695.

**Problem.**
- **Runtime is never given.** At an optimistic 10⁶ Toffolis/s, 10¹¹ Toffolis is about 1 day, 1.3×10¹⁴ (D1
  I2-BDP) is about 4 years, and 6×10¹⁴ (D1 I2-BDP-0) is about 19 years per dye. A reader needs these numbers to
  judge "screening".
- **A single circuit is very deep.** To reach 20% error MLAE needs a deepest Grover power of about 0.6/θ. That
  is m = 2048 for I2-BDP (θ = 3.5×10⁻⁴), or 8×10¹¹ Toffolis in one coherent circuit, and m = 8192 for I2-BDP-0,
  or 3×10¹² Toffolis. The logical error per Toffoli then has to be below about 10⁻¹³. This depth is set by
  1/√p, not by the "2×10⁸ per call" quoted in the abstract.
- **158 qubits is the walk only.** The figure in `tab:decision` covers the THC walk alone. It excludes:
  - the sum-of-Slaters preparation register (log₂D = 14 plus its QROM output);
  - the H_SO block-encoding ancillas;
  - the QSVT signal and LCU qubits.
- **"Required in every case" is too strong.** For the light-atom negatives SPRT is only 2.5–5× more expensive
  than MLAE: BDP 4.6×10¹¹ vs 9.8×10¹⁰, aza-BDP-Ph₄ 1.3×10¹¹ vs 5.2×10¹⁰, BDP-0 1.5×10¹¹ vs 5.6×10¹⁰. In
  addition, SPRT is credited with the true p.

**Fix.**
- Add a runtime column, or a sentence with a stated Toffoli rate.
- Give the maximum coherent depth per circuit and the logical error rate it requires.
- Give the total logical qubits for one call.
- Replace "amplitude estimation is required in every case" with "amplitude estimation is cheaper in every
  case, by 2.5–5× for light-atom negatives and by 10²–10⁴ for rate estimates".

### Q9. The comparisons with Zhou et al. and Eklund et al. need to be tightened

**Where:** l. 80–86, 541–549, 595–600; Fig. `fig:decision`; abstract.

**Problems.**
1. **The Zhou numbers cannot be checked from the repository.**
   - Intro "10⁷–10⁹ Toffoli per shot for 11–45 orbitals" (l. 84–85). The comparison notes
     (`analysis/xanadu_comparison.md`) give 2.7×10⁷ up to 5.8×10⁹ for N = 45, so the upper end is about 10¹⁰.
   - The "483-shot" ISC estimate and the 7.96×10⁷ per shot are hard-coded in `decision_cost.py` (`XANADU`)
     with no table or equation reference.
   - The notes were taken from arXiv v1, while the paper cites the PRA version.
2. **"Like for like" (sec:resources title) is not quite like for like.**
   - **Molecule and active space.** Zhou's N = 11 case is Orig-BDP in the Momeni–Brown π space with frozen
     PCM, and its λ comes from BLISS-THC. Ours is the BODIPY core in the canonical window, which the paper
     itself shows is not a π space, with a double-factorization λ and an unvalidated rank.
   - **Per-shot comparison.** Zhou's 8×10⁷ per shot is essentially their S1 preparation (their H_SO evolution
     is fast-forwarded). It should be compared with our full call (1.2–2.1×10⁸) or with our C_S1 (6.6×10⁷), not
     with the filter alone (l. 547–549).
3. **The abstract omits the proxy comparison [open from previous report, M9].** A reader will compare
   "10¹¹–10¹⁴" with Zhou's 3.8×10¹⁰ without the caveat that the proxy reaches no decision.
4. **The Eklund comparison omits two differences.** "Overlap the lower end" (l. 599–600) is now right
   numerically (Eklund: 10¹³–10¹⁶ Toffolis on 3000–8000 logical qubits). But Eklund et al. simulate the
   nuclear-electronic dynamics, while here the vibronic line shape is a classical MLJ model with assumed
   parameters. That is why our numbers are lower, and it should be said. Our K′ is also exactly the observable
   norm λ_O of Eklund's O(λ_O/ε) accounting; making that link explicit would strengthen the paper.
5. **"Absolute, mechanism-resolved rate" (l. 598) overstates.** The rate is absolute only up to the MLJ
   parameters (λ_s, ħω, S, set a priori, l. 164–165). It also uses vertical gaps at the S0 geometry inside an
   MLJ expression written for adiabatic gaps. The quantum computer improves the electronic part only.

**Fix.**
- Cite table and equation numbers of the PRA version for every Zhou number, and correct "10⁷–10⁹".
- Retitle sec:resources or list the differences.
- Add the proxy caveat to the abstract.
- State the classical line-shape and vertical-gap limitation in sec:filter-method and in the Eklund
  comparison.

### Q10. Cost numbers are tied to an active space and orbitals that the paper itself rejects

**Where:** l. 558–567, `tab:decision` caption, `tab:window`.

**Problem.**
- **The same dye and window have two α_SO values.** `alpha_so.py` uses B3LYP Kohn–Sham canonical orbitals.
  Table `tab:window` reports α_SO "with Hartree–Fock orbitals". For the same dye and window the two differ by
  up to 1.5×:

  | Dye | α_SO, `tab:window` (HF orbitals) | α_SO, `tab:decision` (KS orbitals) |
  |---|---|---|
  | BODIPY core | 133 | 199 |
  | BDP | 27 | 39 |
  | I2-BDP | 8315 | 12 006 |
  | I3-BDP | 11 421 | 12 372 |

  The text says 850–12 400 cm⁻¹ (l. 559) and also "8300 cm⁻¹ in the canonical window" for I2-BDP (l. 564),
  without saying which orbitals each value uses.
- **The active-space dependence is larger still.** It reaches a factor 40 for I2-BDP (299 cm⁻¹ in the AVAS π
  space), and p ∝ α_SO⁻².
- **The costs refer to an uncontrolled window.** The paper argues (l. 646–650) that the right active space
  contains the SO-active orbitals deliberately. The costed one is the uncontrolled canonical window.

**Fix.**
- State the orbitals in the `tab:decision` caption and at l. 559–564, and use one orbital set throughout.
- Present the cost as a function of α_SO/V_eff, with the dyes as illustrative points. Give one variant with an
  AVAS π space plus the iodine 5p (or C–X σ*) orbitals, which is the active space the Discussion recommends.

---

## Minor points

1. **l. 53–57, 691–695: headline ranges.**
   - "Up to 10¹⁴ near the triplet-yield threshold" should read 2–4×10¹⁴, or "a few 10¹⁴".
   - "Between 10¹¹ and 10¹⁴" (Conclusions) does not cover the D1 MLAE range of 5×10¹²–6×10¹⁴ or the D2
     minimum of 5×10¹⁰. Use "5×10¹⁰–6×10¹⁴ (MLAE)" and give the rigorous range (Q5).
   - "10³–10⁶ calls": the source gives 4×10²–3×10⁶. **[open in part from previous report, M6]**
2. **l. 550, `tab:decision` caption: cost per call.** "One call costs 1.4–2.1×10⁸". The D2 calls for clear
   negatives cost 1.19×10⁸, because the polynomial precision is clipped at 0.01, so the range is 1.2–2.1×10⁸.
   (The D1 per-call cost for the light-atom dyes, which are not listed, reaches 3.2×10⁸.)
3. **l. 547, filter walk steps.** "About 3×10⁴ walk steps". The cost model uses 3.2–5.0×10⁴ for D1, which is
   about 1.0–1.5×10⁸ Toffolis for the filter.
4. **l. 544–545, angle-batching overhead.** "For 10–30% more Toffolis". `thc_pennylane_resources.md` gives
   +2% (2 batches), +7% (5) and +17% (10) at optimal swap depth, and up to +25% with swap depth 1. Quote the
   158-qubit point (+7%) explicitly. The rank-44 Toffolis (3.2–5.5×10³) are not given in the text.
5. **l. 511–513, triplet distance.** "Upper triplets that happen to lie within 0.06–0.25 eV". This holds for
   the core and BDP (−0.057, +0.056, +0.065 eV) and for I2-BDP CAS(6,6) (+0.25 eV). For I2-BDP CAS(4,4) the
   nearest upper triplet is at +0.45 eV and the rate comes from Franck–Condon tails. Qualify.
6. **l. 196–204, 570–574: what the simulation tests.**
   - The simulated circuit prepares χ/‖χ‖ exactly (`StatePrep`) and uses a dense unitary dilation of
     (H − E_c)/α with α the spectral range. It does not use an H_SO block encoding, an S1 preparation or an
     LCU/THC encoding of H.
   - "The actual filter circuit" (l. 571) and p = 5.6×10⁻⁶ therefore refer to this normalized circuit. In the
     fault-tolerant model the same dye has p = 1.6×10⁻¹⁵.
   - Say that the simulation verifies the QSVT phases and the polynomial, not the block encodings.
   - Say that the simulation controls the whole QSVT sequence (`qml.ctrl`) while the cost model assumes the
     phase-sign trick with an uncontrolled walk.
7. **Conclusions l. 689–690, "verified to 0.03%".** Add "(2% for the smallest rate)" here as in Results, or
   drop the number. **[open from previous report, M9]**
8. **`tab:filter`, polynomial degree.** The degree of 8000 for I2-BDP CAS(6,6) is the `max_degree` cap of
   `filter_poly`. The error (1.1×10⁻⁵) is still below the tolerance; mention it. A relative-error column was
   requested before and is still missing.
9. **Dynamic-range floor.** The repository also contains an H₂CO run (`H2CO.filtered_isc.json`) where the
   QSVT/polynomial rate (1.2×10⁻⁴ s⁻¹) exceeds the exact rate (4×10⁻¹⁴ s⁻¹) by 10 orders of magnitude. This is
   the dynamic-range floor of about ε²‖χ‖²ρ_max. Either report it as an illustration of the floor or remove the
   file. The ε = r√p/4 rule handles the floor correctly. Derive the rule in the SI:
   δp ≤ 2ε√(qp) + ε²q with q = ‖χ‖²/α_SO² ≤ 1. The bound shows the rule is conservative by 1/√q.
10. **l. 228–231, assumptions.** γ² = 0.5, D = 10⁴ determinants and the 0.5 eV S1–S2 gap are still unjustified.
    Take them from the CAS(12,11) S1 vector (the number of determinants for 50% weight) and the CAS S1–S2 gap.
    **[open from previous report, minor 14]**
11. **`c_so` (H_SO basis rotation).** The function uses N(N−1)/2 Givens rotations for N = 11 spatial orbitals.
    The SOC matrix is a complex 22×22 spin-orbital matrix, so 2N(2N−1)/2 complex rotations plus the
    Σμ_k Z_k LCU are needed. The cost is still negligible, but state it correctly.
12. **QSVT phase rotations.** The QSVT phase rotations (one Rz per step at 27 bits) and the reflections are
    not in C_call. They are about 1% of the walk cost; say so.
13. **QPE confidence.** The QPE cost uses πλ′/(2ε) walk steps, which is an RMS-type precision. A 95% bound
    needs a constant factor more. This is negligible in the totals; state the confidence.
14. **`tab:decision` scope.** The table lists 8 of 15 dyes, but the text ranges ("light-atom negatives
    5×10¹⁰–5×10¹¹") come mostly from dyes that are not in it (BODIPY core, aza-BDP-Ph₄, BDP-0, Cl2-BDP,
    azaBDP_BrPh2). Put the full table in the SI; this is in the SI TODO at l. 709.
15. **SI Fig. S3 percentiles.** The crosses (`circuit_qpe`) are 90th-percentile absolute errors divided by p.
    The curves are 95th-percentile relative errors. Use the same percentile.
16. **l. 578, wording.** "Neither amplitude estimator resolves p below about π/√p ≈ 1.3×10³ Grover
    iterations" should read "… needs at least about π/√p … Grover iterations before its error falls below
    100%".
17. **`decision_cost.py`, D2 MLAE confidence.** `mlae_calls` for D2 uses the two-sided 95th percentile as a
    one-sided bound, so the per-decision confidence is about 97.5%. This is slightly conservative; state it.

---

## Numbers checked

The sources below are under `paper/quantum/` or `psscreen/`. Numbers not listed under Major or Minor points
match their source.

- **Cost-model constants** (`decision_cost.json` model).
  - C_walk 2.97×10³ and 158 qubits.
  - λ′ = 15.9 + 3.46 = 19.4 Ha, from `degree_check.json`: |E0 − c0| = 3.4607 Ha.
  - Δ = √(2λ_s k_BT) = 0.0717 eV.
  - Degree fit 0.554 ln(1/ε) − 0.106.
  - C_S1 6.58×10⁷, C_SO 4.5×10³, C_QPE 1.11×10⁸.
  - N_rep 16 = ⌈ln 20 / (2(8/π² − ½)²)⌉.
  - c_MLAE 45.04, the median of c_eff for K = 7–13.
- **eq `eq:p`, D1.**
  - I2-BDP: ε = r√p/4 = 1.7×10⁻⁵ gives degree 43 779 (checked by hand).
  - DKLR Y₁ = 319 gives 2.6×10⁹ calls.
  - MLAE 45/(0.2 × 3.46×10⁻⁴) = 6.5×10⁵ calls.
  - Total 1.3×10¹⁴.
- **Table `tab:decision`.** All 8 rows (α_SO, ‖χ‖, k/k_th, p, D1 MLAE, D1 sampling, D2 best) match
  `decision_cost.md`.
- **sec:decision.**
  - α_SO 850–12 400 cm⁻¹ (KS orbitals).
  - α_SO/‖χ‖ = 5.2–156 for the ISC-active dyes; 728 and 668 for I-BDP and I-BDP-0.
  - p = 6.2×10⁻⁹–4.4×10⁻⁵.
  - D1 MLAE 5.4×10¹²–6.1×10¹⁴; ratio to one call 3.4×10⁴–2.9×10⁶; sampling/MLAE 213–17 977.
  - D2: light-atom negatives 5.2×10¹⁰–5.4×10¹¹; heavy-atom dyes far from the threshold 5.8×10¹¹–1.7×10¹³;
    Br2-BDP 2.3×10¹⁴ and azaBDP_IPh2 3.7×10¹⁴, at 0.94× and 1.2× the threshold.
  - Family-wise factor ln 4000/ln 40 = 2.25.
  - Zhou total 7.96×10⁷ × 483 = 3.8×10¹⁰; D1 ratio 1.4×10²–1.6×10⁴.
- **AE demo** (`ae_demo.json`).
  - Grover law to 5.4×10⁻¹⁰.
  - p = 5.56×10⁻⁶ (M_S = 1 sector).
  - MLAE: 17.9% at 1.03×10⁵ calls, 9.3% at 2.06×10⁵.
  - Canonical: 17.6% at 2.62×10⁵ (ensemble).
  - Sampling: 248% and 162%.
  - π/√p = 1332; c_canonical 109.
- **Filter simulation** (`*.filtered_isc.json`, `tab:filter`).
  - Every QSVT, exact, V1 and ΔE entry matches.
  - QSVT/exact agree to 0.03% for I2-BDP (6,6) and 2.4% for I2-BDP (4,4).
  - Wires: 6/8 and 10/11.
  - Degrees 3324–8000.
  - At CAS(6,6) the ranking is reversed: core 2.48×10⁸ vs I2-BDP 4.23×10⁶ s⁻¹.
- **THC** (`thc_pennylane_resources.md`).
  - Coefficient bits 20 and rotation bits 27 (from the stated rules with ε = 0.16 mHa).
  - Rank 22: 2.77–3.46×10³ Toffolis per step; 374/239/158/131 wires.
  - Rank 44: 136–379 wires.
- **Re-run** (`degree_check.min_degree`). Degrees 11 028 (λ/Δ = 1600, ε = 10⁻⁵), 11 868 (1600, 3×10⁻⁶) and
  2756 (400, 10⁻⁵), against fit values of 10 035, 11 103 and 2509.
- **Derived for this report** (from `decision_cost.json`).
  - V_eff = √p·α_SO: 4.2 cm⁻¹ (I2-BDP), 5.7 (azaBDP_Br2), 0.21 (I2-BDP-0), 1.0 (azaBDP_IPh2).
  - Canonical-QAE totals 6–8× MLAE.
  - Coherent MLAE depth at 20%: m = 2048 (I2-BDP) and 8192 (I2-BDP-0).
