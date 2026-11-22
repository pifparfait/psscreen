# Decision-aware resource estimate (quantum contribution, option A), revised after review

> **Superseded (2026-11-22).** These are working notes from the first version of the cost model. The current model (λ = 17.3 Ha regularized THC, filter degree computed directly by `degree_direct.py`, QPE at 0.1 mHa, D2 with a factor-2 indifference zone and no known margin) is documented in the docstring of `decision_cost.py`, and its numbers are in `results/decision_cost.md`.


Code: `paper/quantum/decision_cost.py` (cost model and tables), `alpha_so.py` (spin-orbit block-encoding
normalization), `degree_check.py` (filter degree from our own polynomials), `ae_demo.py` (amplitude estimation on the
actual filter circuit). Results: `paper/quantum/results/decision_cost.md`, `alpha_so.json`, `degree_check.json`,
`ae_demo.json`; figures `paper/figures/decision_cost.pdf`, `ae_demo.pdf`. Review that triggered this revision:
`paper/analysis/decision_cost_review.md`.

## Question

Previous resource estimates for photosensitizer screening (Zhou et al., Phys. Rev. A 2026, arXiv:2512.15889) report
Toffolis per shot, plus shot counts for a fixed additive precision of a relative proxy. A screening campaign needs a
decision per molecule: "is Φ_T above 0.5?" or "k_ISC within the factor-2 budget". What does that decision cost?

## Model (what changed in the revision is marked ★)

- ★ **Success probability with the true normalization.** χ = H_SO|S1⟩ is produced by a block encoding of the
  one-body H_SO, whose normalization is α_SO = ½ Σ|μ_k| (eigenvalues μ_k of the spin-orbital SOC matrix; one basis
  rotation). One call succeeds with p = k_ISC/K′, K′ = (2π/ħ) α_SO² ρ_max. K′ is known exactly, so the previous
  ‖χ‖²-from-6-triplets caveat is gone.
- Cost per call, CAS(12,11), THC qubitization (158 logical qubits, 2.97×10³ Toffolis per walk step):
  ★ coherent S1 preparation (sum of Slaters with 10⁴ determinants, an S1 eigenstate filter for a 0.5 eV S1–S2 gap,
  one round of amplitude amplification for an assumed overlap γ² = 0.5; 6.6×10⁷), the H_SO basis rotation, and the
  filter. ★ Filter degree from our own Chebyshev filters: d = (0.55 ln(1/ε) − 0.11) λ′/Δ, with
  ★ λ′ = λ + |E₀ − c₀| = 15.9 + 3.5 = 19.4 Ha (shift relative to the identity coefficient). The preprint-derived
  constant used before overestimated the degree by 1.3–2×; the shift adds 22%. Net: 1.6–2.1×10⁸ Toffolis per call
  for the ISC-active dyes. One QPE for E(S1) at 2 mHa (1.1×10⁸) is added once; its error is a systematic gap error.
- ★ **Estimators on equal terms**, confidence 0.95 per decision:
  sampling with a stopping rule that does not need p (Dagum–Karp–Luby–Ross) for D1 and Wald's SPRT for D2;
  maximum-likelihood amplitude estimation (MLAE: no QPE, no controlled walk) with the empirical constant from the
  simulation; canonical QAE (power-of-two M, median of 16) as a rigorous check.

## Findings

1. ★ **The spin-orbit normalization is the dominant cost driver.** For heavy-atom dyes α_SO is 850–12 400 cm⁻¹
   in CAS(12,11), 7–150 times ‖H_SO S1‖ (TDDFT). The success probability of the ISC-active dyes is therefore
   p = 6×10⁻⁹ to 4×10⁻⁵, 10²–10⁴ times smaller than in the first version. Most of the norm of H_SO couples S1 to
   states far from resonance, or is never touched by S1.
2. **Per-call costs understate the decision cost by 10⁴–10⁶.** k_ISC to 20% costs 5×10¹²–6×10¹⁴ Toffolis with
   MLAE (3×10⁴–3×10⁶ calls of 1.6–2.1×10⁸). The ratio total/per-call is 3×10⁴–3×10⁶ for D1 and 4×10²–2×10⁶ for D2.
3. ★ **Amplitude estimation is required for every decision.** Direct sampling costs 2×10²–2×10⁴ times more for D1
   (10¹⁵–10¹⁹ Toffolis). For D2, MLAE is also cheaper for every dye, including clear negatives, because certifying
   p < p_th takes about 1/p_th samples but only about 1/√p_th amplitude-estimation calls, and p_th is small once α_SO
   is in the denominator. The "sample clear negatives" part of the earlier adaptive strategy does not survive.
4. **The cost of a classification is set by the margin to the threshold.** Light-atom negatives: 5×10¹⁰–4×10¹¹
   Toffolis. Heavy-atom dyes far from Φ_T = 0.5 (both sides): 6×10¹¹–2×10¹³. Dyes near the threshold
   (Br2-BDP at 0.94× the threshold rate, azaBDP_IPh2 at 1.2×): 2–4×10¹⁴.
5. ★ **Simulation on the actual filter circuit** (BODIPY core, CAS(4,4), M_S = 1 sector, p = 5.6×10⁻⁶): the
   good-state probability after m Grover iterations follows sin²((2m+1)θ) to 4×10⁻⁹, so QAE and MLAE can be
   simulated from θ. At the 95th percentile: MLAE reaches 18% relative error with 1.0×10⁵ calls and 9% with 2.1×10⁵;
   canonical QAE (median of 16, averaged over nearby amplitudes) needs 2.6×10⁵ for 18%; direct sampling is still at
   250% and 160%. The empirical constants (ε_rel·√p·N) are 45 for MLAE and 109 for canonical QAE. The earlier
   "5% at 1.3×10⁵ calls" was a grid coincidence (2ᵗθ/π = 3.07) and is withdrawn. Neither estimator resolves p below
   about π/√p ≈ 1.3×10³ Grover iterations. With 100 shots per depth MLAE flattens at about 3% (rare branch errors),
   so the constant is used only for r ≥ 0.05.
6. **Lever for future work.** Because the cost scales as α_SO/‖χ‖ (amplitude estimation), a tighter H_SO encoding
   pays off directly: for example, encoding only the heavy-atom one-centre part, or preparing χ by a filtered
   (resonance-weighted) state preparation. Quantifying this is out of scope.

## Paper phrasing (replaces "first decision-level resource estimate")

"Previous resource estimates for photosensitizer screening report costs per shot, or per fixed additive precision of
a relative proxy [Zhou2026PRA]. Following end-to-end accountings that include the measurement stage
[KnillOrtizSomma2007, Steudtner2023, Eklund2026], we estimate the cost of the screening decision itself. The
success probability of the golden-rule circuit is k_ISC/K′ with K′ fixed by the spin-orbit block-encoding norm;
for heavy-atom dyes it is 10⁻⁸–10⁻⁵, so the number of circuit calls, not the cost per call, dominates. Amplitude
estimation is required for every decision, and the cost of a classification is set by how close the dye lies to the
threshold: 10¹¹ Toffolis for clear negatives, 10¹³ for clear positives and 10¹⁴ near Φ_T = 0.5."

Comparison with the preprint: our decisions cost 10²–10⁴ times the preprint's 3.8×10¹⁰ Toffolis for the ISC proxy
(483 shots). This is not like for like. The proxy targets a fixed additive precision of a relative, gap-free quantity
and cannot take decision D2 without calibration; the difference is the price of an absolute, mechanism-resolved rate.
Our totals fall in the range Eklund et al. report for photochemical yields (10¹³–10¹⁶ Toffolis).

## Assumptions to state in the paper

- Electronic inputs (k_ISC, k_r, k_ic) are the TDDFT/ΔSCF pipeline values; α_SO is from the CAS(12,11) SOC
  integrals of the same dye.
- S1 overlap γ² = 0.5 and a 0.5 eV S1–S2 gap for the eigenstate filter are assumptions.
- r = 0.2 uses about a quarter of the factor-2 budget in ln k, leaving the rest to the systematic gap error;
  r = 0.4 halves the D1 costs.
- δ = 0.05 per decision; for a family-wise error over N dyes use δ/N (×2.2 for N = 100).
- Canonical QAE omits the cost of controlling the walk operator; MLAE needs no control.
