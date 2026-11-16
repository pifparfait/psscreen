# Review of the decision-aware resource estimate (commit 6fa9947)

Scope: `paper/quantum/decision_cost_notes.md`, `decision_cost.py`, `ae_demo.py`, `results/decision_cost.md` and
`results/ae_demo.json`. This is a novelty check against the literature plus a referee-style critique. No code was
changed. New references are appended to `paper/related_work.bib` under "DECISION-LEVEL COST AND AMPLITUDE
ESTIMATION". Each one was checked against its Crossref DOI record or its arXiv listing on 2026-11-16.

## Summary verdict

The idea is sound and worth keeping: the observable's success probability is intrinsically small, so the call count
dominates the per-call cost, and the cost of a classification is set by the margin to the threshold. Four issues
change numbers or claims and should be fixed before the paper uses them:

1. **"First" cannot be claimed as written.** End-to-end costs that include amplitude estimation for photochemical
   observables (quantum yields and rates) were published in March 2026 by Eklund et al. (Google/Sydney). Observable
   estimation with amplitude estimation is standard in fault-tolerant chemistry costing.
2. **The success probability omits the block-encoding normalization of H_SO.** p is defined relative to ‖χ‖², but
   the circuit cannot prepare χ/‖χ‖ for free. This is the largest correction. It probably lowers p by one to two
   orders of magnitude, and it also removes the ‖χ‖² truncation caveat.
3. **The 5% amplitude-estimation result in the demo is a grid coincidence, not typical performance.**
4. **Some headline numbers are inconsistent with the generated table.** The "10³–10⁵" understatement factor and the
   "π/(2√p)" threshold both disagree with it.

Details and concrete fixes follow, in order of importance.

---

## 1. Novelty: what exists, and what can be claimed

| Work | What it already does | Consequence for our claim |
|---|---|---|
| Eklund, Tikku, Sinnott, Huggins, Low, Berry, Kassal, arXiv:2603.19007 (Mar 2026), `Eklund2026` | End-to-end fault-tolerant algorithm for photochemical dynamics. It covers state preparation, qubitized evolution and **measurement of yields and rates with amplitude estimation** (query cost O(λ_O/ε)). It quotes total costs for a quantum yield to ε = 0.095: 10¹³–10¹⁶ Toffolis on 3000–8000 logical qubits. | Total cost including amplitude estimation for a photochemical rate or yield is **not new**. Cite it, and note that our totals of 10¹²–10¹⁴ Toffolis fall in the same range. |
| Knill, Ortiz, Somma, PRA 75, 012328 (2007), `KnillOrtizSomma2007` | Expectation values with Heisenberg (1/ε) scaling via amplitude or phase estimation. | Origin of "amplitude estimation is not optional"; cite it. |
| Rall, PRA 102, 022408 (2020), `Rall2020` (already in the .bib) | Spectral densities and other physical quantities from block encodings, with amplitude-estimation-based cost. | Same primitive as ours: the filter plus amplitude estimation. |
| Huggins et al., PRL 129, 240501 (2022), `Huggins2022` | Nearly optimal estimation of many expectation values at the same time. | Relevant if several dyes or channels are estimated together. |
| Steudtner et al., Quantum 7, 1164 (2023), `Steudtner2023` (already in the .bib); Apel, Cortes, Lemieux, Steudtner, arXiv:2508.06677 (2025), `Apel2025` | Fault-tolerant cost of molecular observables, with total Toffolis per observable at a target precision. Apel et al. use windowed coherent QPE and report a >100× Toffoli reduction. | Per-observable totals with a 1/ε estimator are standard in PsiQuantum's chemistry resource estimates. |
| Lin & Tong, PRX Quantum 3, 010318 (2022), `LinTong2022` | Ground-state energy estimation cast as a CDF or decision problem, with Heisenberg-limited cost. | Precedent for "cost of a decision set by the distance to a threshold", but for energies. |
| Zhou et al., Phys. Rev. A (2026), `Zhou2026PRA` | Gives shot counts (483 for the ISC proxy, 133 for the window) next to per-shot costs. | Their totals are not "per-shot only"; they are an additive-error total for a different observable. See §2.6. |

**What is new, as far as I found:**

- For a golden-rule spectral-filter observable, the success probability scales as k_ISC / K. That makes the call
  count, not the per-call cost, the dominant term for exactly the dyes that matter.
- For PDT screening, we use a classification criterion (Φ_T > 0.5) whose cost depends on the margin to the
  threshold. This gives an adaptive strategy (sample first, then switch to amplitude estimation).
- The cost is shown to be 2–3 orders of magnitude above the per-observable totals quoted for PDT screening so far.

**Suggested phrasing (replaces "the first decision-level resource estimate"):**

> "Previous resource estimates for photosensitizer screening report costs per shot or per fixed additive precision
> [Zhou2026PRA]. Following end-to-end accountings that include the measurement stage [KnillOrtizSomma2007,
> Steudtner2023, Eklund2026], we estimate the cost of the screening decision itself. Because the golden-rule rate is a
> small fraction of the total spin–orbit weight, the number of circuit calls dominates. The cost of a classification
> is set by how close the dye lies to the decision threshold."

---

## 2. Errors and weak points in the model

### 2.1 Missing H_SO normalization in p (major)

`isc_filtered.py` and `ae_demo.py` prepare χ/‖χ‖ with `StatePrep` (`phi[:dim] = c`, with `c` normalized). On a
fault-tolerant machine, χ = H_SO|S1⟩ has to be produced by a **block encoding of H_SO**. Its normalization is α_SO,
for example the sum of absolute eigenvalues of the three one-body SOC components, or the LCU one-norm. One call of A
then succeeds with probability

  p_tot = ⟨S1|H_SO P(H)² H_SO|S1⟩ / α_SO² = p · ‖χ‖²/α_SO²,

and the rate is

  k_ISC = K′ p_tot,  with K′ = (2π/ħ) α_SO² ρ_max.

The basis-rotation cost `c_so()` covers the unitary part. It does not change this success probability, which is
what the call count depends on. Consequences:

- **Costs go up.** The amplitude-estimation cost scales by α_SO/‖χ‖ and the sampling cost by (α_SO/‖χ‖)².
  Amplifying the χ preparation instead gives the same α_SO/‖χ‖ factor, applied to the per-call cost. I could not
  compute α_SO here (no PySCF in this container). With iodine p orbitals in the active space, α_SO is plausibly
  thousands of cm⁻¹. Our ‖χ‖ is about 10³ cm⁻¹ for I2-BDP (from K = 4.9×10¹⁴ s⁻¹) and about 3 cm⁻¹ for the BODIPY
  core. So ‖χ‖²/α_SO² ≈ 10⁻²–10⁻¹ is a reasonable guess, but **it must be computed**.
- **The ‖χ‖² truncation caveat disappears.** With K′, the circuit needs α_SO, which is known exactly from the
  Hamiltonian, not ‖χ‖², which was estimated from 6 TDDFT triplets. The ×10 sensitivity column then becomes
  unnecessary.
- **Some negatives cost nothing.** If K′ < k_th, then k_ISC ≤ K′ < k_th and the classification is decided without
  any circuit call. The current table already has a case like this: BDP-0 has p_th = 1.28 > 1, so the reported
  2×10⁸ Toffolis should be 0. Under K′ fewer dyes will qualify, but the check is worth adding.

**Fix:** compute α_SO for the CAS(12,11) SOC operator, define p_tot and K′ as above, and rerun both tables. In
`ae_demo.py`, either prepend a block encoding of H_SO or multiply p by ‖χ‖²/α_SO² and say so.

### 2.2 The amplitude-estimation demo: the 5% result is a grid coincidence (major for the text, not for the model)

In `ae_demo.json`, at t = 12 and t = 13 the AE error is exactly 4.8%, both at the median and at the 90th percentile.
Canonical QAE returns p̂ = sin²(πy/2ᵗ). For p = 5.56×10⁻⁶ we have θ = 2.36×10⁻³ and 2ᵗθ/π = 3.07 at t = 12, so the
nearest grid point y = 3 happens to sit 4.8% from p. The median of 16 runs does not remove this bias, because every
run lands on the same grid point. For a generic θ at the same t, the error can be anywhere up to about 2·(0.5/3), or
33%.

**Fix:**
- Average the demo over a range of p, for example by shifting the filter centre slightly or scaling the polynomial
  amplitude, or report the worst case over grid offsets.
- Better, simulate a grid-free estimator (MLAE or IQAE, §3). The current text "reaches 5% relative error with
  1.3×10⁵ calls" overstates canonical QAE.
- The model itself uses the Brassard bound, which already covers grid error, so `decision_cost.py` is not affected.

Also, AE first resolves p when 2ᵗ ≳ π/θ = π/√p (y ≥ 1), not at π/(2√p). With p = 5.6×10⁻⁶ that is about 1.3×10³
iterations, which is consistent with the t = 11 → 12 transition in the data. Please correct that sentence.

### 2.3 Headline numbers that disagree with the generated table

- "Per-shot counts understate the cost of a decision by 10³–10⁵." From `decision_cost.md`, D1 with AE costs
  2.4×10¹²–6.4×10¹³ Toffolis against 1.0–1.2×10⁸ per call, a ratio of 2×10⁴–5×10⁵. With sampling the ratio is
  10⁵–5×10⁷. D2 runs from about 2× for clear negatives to 10⁵ near threshold. State the ratio per decision type,
  or say "10⁴–10⁵ for rate estimates with amplitude estimation".
- The column "Grover calls (AE)" actually lists `ae_rep*(2M+1)`, which is calls to A or A†. Rename it.

### 2.4 Estimator bounds: amplitude estimation and sampling are not compared on equal terms

- **Canonical QAE (Brassard et al. Thm 12).** The bound |p̂ − p| ≤ 2π√(p(1−p))/M + π²/M² with probability ≥ 8/π² is
  correct. The median boost n_rep = ⌈ln(1/δ)/(2(8/π² − ½)²)⌉ = 16 for δ = 0.05 is correct, by Hoeffding on the
  majority. Two details:
  - Canonical QAE needs M = 2ᵗ. `ae_iterations` returns any integer M, which can underestimate by up to 2×. Round
    up to a power of two, or say that a grid-free estimator is assumed.
  - Calls per run should be 2(M−1)+1 controlled operations, not 2M+1. That is minor, but each Grover iterate is
    *controlled*, and controlling the THC walk adds per-step overhead. That overhead is not in C_call. Either add
    it or note that it is neglected.
- **Bernstein sampling.** n = 2p(1−p)L/ε² + 2L/(3ε) is the correct inversion of Bernstein's inequality. However,
  it uses the **true** p. A real relative-error procedure does not know p and needs a sequential stopping rule. The
  Dagum–Karp–Luby–Ross rule (`Dagum2000`) costs about 1.5–2× more. Combined with the 16× median and the 2π
  constant on the AE side, both columns sit at worst-case constants of different looseness. The quoted 4–92×
  advantage is therefore **probably an underestimate** of the practical AE advantage (§3).
- **D2 should be a sequential test, not a fixed-precision estimate.** For a one-sided decision against p_th, a Wald
  SPRT (`Wald1945`) needs about ln(1/δ)/KL(p‖p_th) samples. For I2-BDP (p/p_th ≈ 20) that is about
  3/(1.4×10⁻⁵ × 2.0) ≈ 1×10⁵ shots, against the 7.6×10⁵ in the table. Clear positives are therefore over-costed
  by about 5–10× on the sampling side. Near-threshold dyes are unaffected, because there the margin dominates.
- **Campaign-level confidence.** δ = 0.05 *per decision* gives about one misclassification per 20 dyes. For N dyes,
  use δ/N (union bound). The cost grows only logarithmically, but the choice should be stated.

### 2.5 Degree model

- Using the preprint's step-function fit for a Gaussian-type filter is a reasonable heuristic. Gaussian filters
  scale as (λ/σ)√ln(1/ε), so the √ln form is right. Our own simulated degree is a quick sanity check: CAS(4,4)
  BODIPY core gives d = 3324 at tol = 2×10⁻³. Fit d against α/Δ from `filter_poly` over several α values and
  report that fit next to the preprint-derived one. If they agree within about 20%, say so. That turns a borrowed
  constant into a validated one.
- **λ versus λ′.** The preprint's fit uses λ′ = λ + |E_th|. Here (H − E_c)/α is block-encoded with the THC λ =
  15.9 Ha, but the shift by E_c adds |E_c| (in the block-encoding frame) to the normalization. If |E_c| is
  comparable to λ, the degree and the per-call cost roughly double. Check the frame and use λ′ = λ + |E_c|.
- The choice ε_poly = r√p/4 follows from |p̃ − p| ≤ 2ε√p + ε². That is correct. Note that after §2.1, p should
  become p_tot, which makes ε_poly smaller. The effect is only on the √ln factor.

### 2.6 S1 preparation and the E(S1) QPE

- **Exact S1 is assumed in every call.** Amplitude estimation needs A to be a fixed unitary, so S1 preparation must
  be coherent in every call. Any S0, S2, … admixture of amplitude c contributes cross terms of order
  2|c|·V_S0,T·V_S1,T·ρ to p. Because p is tiny, a 1% admixture can bias it by more than the 20% target if the
  contaminant couples strongly to resonant triplets. Two options:
  - (a) Add a coherent S1 eigenstate filter to each call: degree about (λ/ΔE_S2−S1)·ln(1/ε), roughly 3×10³ walk
    steps for a 0.8 eV S1–S2 gap. That adds about 10% to C_call and multiplies p by the overlap γ².
  - (b) State the overlap requirement γ² ≥ 1 − O(r·p_res) explicitly.
  Either way, the "sum of Slaters with 10⁴ determinants, assumed exact" caveat should move from the caveats into
  the model.
- **E(S1) QPE once at 2 mHa.** Using π λ/(2ε) walk steps is fine. This QPE fixes the filter centre, so its error
  enters as a systematic gap error, not as statistics. That is consistent with the accuracy budget (2 mHa per state
  is about 0.05 eV). Say so, and note that the QPE also needs a good S1 overlap.

### 2.7 Comparison with the preprint (483 shots × 8×10⁷)

That comparison is not like-for-like:
- The proxy targets a fixed **additive** resolution (10%) on an O(1)-scaled probability whose magnitude is set by
  the chosen evolution time.
- The proxy gives no absolute rate. It cannot take decision D2 at all without calibration against a reference
  molecule, which means at least a second estimate with its own error.

Fair phrasing: "the gap-resolved observable costs 10²–10³ times more per decision than the preprint's proxy costs
per observable, because the information it carries (an absolute, mechanism-resolved rate) sits in a small
amplitude." Keep the horizontal line in the figure, but label it "proxy, fixed additive precision (not a decision)".

### 2.8 Smaller points

- **r = 0.2 is stricter than the factor-2 budget.** A factor of 2 corresponds to |ln k| ≤ 0.69, so a statistical
  r ≈ 0.3–0.5 already fits inside the budget. AE cost ∝ 1/r, so using r tied to the budget cuts D1 by 1.5–2.5×. If
  0.2 is kept, justify it as leaving room for systematic error.
- **K varies from 3.5×10⁸ to 8×10¹⁴ s⁻¹ across dyes** only because ‖χ‖² comes from whichever 6 triplets TDDFT
  returned. This is another reason to switch to α_SO (§2.1).
- **D1 restricted to k_ISC ≥ 10⁶ s⁻¹.** Fine, but give the reason in the paper: below that, Φ_T is far from 0.5, so
  only D2 matters.

---

## 3. Better estimators: would they change the numbers?

Yes, in the direction that strengthens the claim.

| Estimator | Key property | Expected effect on our numbers |
|---|---|---|
| Canonical QAE, median of 16 (current) | Power-of-two grid, worst-case 2π constant, 16 repetitions | Baseline. Calls ≈ 32M with ε ≈ 2π√p/M, i.e. ε ≈ 200√p/N_calls |
| MLAE (`Suzuki2020`) | No QPE; Grover depths 0, 1, 2, 4, …, maximum likelihood; near Heisenberg-limited Fisher information | Removes the grid bias of §2.2 and the 16× median. Constant-factor saving, likely 5–20× (simulate to confirm) |
| IQAE (`Grinko2021`) | Rigorous confidence intervals, no QPE, adaptive depth | Same order as MLAE, with a provable δ |
| Approximate counting, simplified (`AaronsonRall2020`) | Simple Heisenberg-scaling estimator with a rigorous bound | Alternative rigorous baseline |
| Low-depth QAE (`GiurgicaTiron2022`) | Trades depth against total calls: queries ∝ ε^{−1−β} at depth ∝ ε^{−1+β} | Relevant because one call is already ~10⁸ Toffolis. Gives a depth-limited curve for early fault-tolerant hardware |
| AE from QSP (`RallFuller2023`) | Near-optimal constants, coherent | Best rigorous constants |
| Relative-error AE (`Montanaro2015`) | O(1/(r√p)) without knowing p in advance | Justifies the relative-error D1 scaling without the true p |
| Windowed coherent QPE (`Apel2025`) | Kaiser windows cut Toffolis >100× for observables | Applies if A is reformulated as an observable estimate |

**Recommendation:** rerun the cost model with a grid-free estimator (MLAE or IQAE) and an SPRT for D2. In
`ae_demo.py`, add an MLAE run on the same circuit and report the empirical constant c in ε ≈ c√p/N_calls. Expect
the AE/sampling advantage to rise from 4–92× to roughly 20–500×. The qualitative conclusions ("AE is not
optional", "cost set by the margin", "adaptive strategy") all survive, and the first one becomes stronger.

---

## 4. Checklist of fixes (in order)

1. Add α_SO normalization: p_tot = p‖χ‖²/α_SO² and K′ = (2π/ħ)α_SO²ρ_max. Drop the ‖χ‖²×10 column. Decide the
   free negatives with K′ < k_th.
2. Fix the demo: average over grid offsets or use MLAE. Correct the "5% at 1.3×10⁵ calls" and "π/(2√p)" sentences.
3. Replace "first decision-level resource estimate" with the phrasing in §1. Cite Eklund2026, KnillOrtizSomma2007,
   Rall2020, Steudtner2023 and Apel2025.
4. Use a sequential test for D2 (SPRT) and a stopping rule for D1 sampling (Dagum2000). Round canonical-QAE M up to
   powers of two, or switch to MLAE or IQAE.
5. Make S1 preparation part of the model: add a per-call eigenstate filter, or state the overlap requirement.
6. Check λ′ = λ + |E_c| in the degree. Validate the degree constant against our own `filter_poly` degrees.
7. Report the understatement factor per decision type. Rename the "Grover calls" column. Justify r = 0.2 or tie it to
   the budget. State δ per decision versus per campaign.
8. Relabel the preprint reference line as "fixed additive precision, not a decision".
