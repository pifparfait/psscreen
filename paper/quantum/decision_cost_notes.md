# Decision-aware resource estimate (quantum contribution, option A)

Code: `paper/quantum/decision_cost.py` (cost model and tables), `paper/quantum/ae_demo.py` (amplitude estimation
simulated in PennyLane). Results: `paper/quantum/results/decision_cost.md`, `ae_demo.json`; figures
`paper/figures/decision_cost.pdf`, `ae_demo.pdf`.

## Question

Published resource estimates, including arXiv:2512.15889 (Phys. Rev. A), give Toffolis **per shot**. A screening
campaign needs a **decision** per molecule, for example "is Φ_T above 0.5?" or "k_ISC to within the factor-2 budget".
How many calls to the circuit does that decision need, and what does it cost in total?

## Model

- The gap-resolved filter returns k_ISC = K·p, where p is the success probability of one circuit A and
  K = (2π/ħ)‖H_SO|S1⟩‖² ρ_max.
- Cost of one call for CAS(12,11) with THC qubitization (158 logical qubits, 2.97×10³ Toffolis per walk step):
  S1 preparation, plus the one-body H_SO by basis rotation, plus a filter of degree d ≈ 3.3-4.1×10⁴. That gives
  1.0-1.2×10⁸ Toffolis, the same order as the preprint's ISC proxy per shot (8×10⁷). One QPE for E(S1) at 2 mHa
  adds 3.7×10⁷ once.
- Number of calls, comparing direct sampling of the success flag (as in the preprint) with canonical amplitude
  estimation (median of 16 runs, confidence 0.95). Sampling scales as 1/(r²p) for relative error r; amplitude
  estimation scales as 1/(r√p).

## Findings

1. **p is small for every dye that matters, so per-shot costs understate the cost of a decision by 10³-10⁵.**
   The ISC rate is a tiny fraction of the total spin-orbit weight: most of ‖H_SO S1‖² goes to triplets far from
   resonance. For the ISC-active dyes p = 4×10⁻⁶ to 2×10⁻³. Estimating k_ISC to 20% then costs
   2.4×10¹²-6.4×10¹³ Toffolis with amplitude estimation, and 1×10¹³-6×10¹⁵ with direct sampling.
2. **Amplitude estimation is not optional.** It reduces the cost of the rate estimate by 4-92× for the 9 dyes
   with k_ISC ≥ 10⁶ s⁻¹. The saving is largest for the iodinated dyes (50-92×), which are the most relevant.
3. **The cost of a decision is set by the distance to the decision boundary.** Clear negatives (no heavy atom,
   Φ_T ≈ 0) are classified with a few to a few hundred shots, 2×10⁸-2×10¹⁰ Toffolis, and direct sampling is
   cheaper there. Dyes near Φ_T = 0.5 (Br2-BDP at 0.94× the threshold rate, azaBDP_IPh2 at 1.2×) need
   6×10¹²-1.5×10¹³ Toffolis even with amplitude estimation (10-23× cheaper than sampling). The best strategy is
   adaptive: sample first, then switch to amplitude estimation when the margin is small.
4. **Verified in simulation.** In PennyLane, canonical QAE on the actual filter circuit (BODIPY core, CAS(4,4),
   M_S = 1 sector, p = 5.6×10⁻⁶) reaches 5% relative error (median and 90th percentile, median of 16 runs) with
   1.3×10⁵ calls. With 2.6×10⁵ calls, direct sampling still has a 37% median error and 106% at the 90th
   percentile. Below about π/(2√p) Grover iterations neither estimator resolves p, which is the expected
   Heisenberg threshold.
5. **Sensitivity.** ‖H_SO S1‖² is computed from 6 triplets per dye, so it is a lower bound. Multiplying it by 10
   raises the amplitude-estimation costs by about 3.4× (√10), and the conclusions do not change.

## How to phrase it in the paper (quantum contribution)

"We provide the first decision-level resource estimate for photosensitizer screening on a fault-tolerant quantum
computer. Because the ISC rate is a small fraction (10⁻⁶-10⁻³) of the total spin-orbit weight, the number of circuit
calls, not the cost per call, dominates; per-shot estimates understate the cost of a screening decision by three
to five orders of magnitude. Amplitude estimation reduces this by up to two orders of magnitude for the iodinated
dyes, and an adaptive sampling/amplitude-estimation strategy keeps clear negatives cheap."

Caveats: the walk-step and degree models are taken from PennyLane's THC estimator and from the preprint's empirical
degree fit (with a Gaussian-filter scaling in the polynomial precision). The electronic inputs are TDDFT/ΔSCF
quantities, and S1 preparation is assumed exact (sum-of-Slaters with 10⁴ determinants, as in the preprint).
