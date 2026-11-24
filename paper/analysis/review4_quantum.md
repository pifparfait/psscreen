# Referee report (quantum algorithms and resource estimation), round 4

Manuscript: "How Accurate Must Excited States Be to Screen Triplet Photosensitizers?" (main.tex, si.tex, si_tables.tex),
cost model in paper/quantum/. Line numbers refer to main.tex and si.tex as of 2026-11-24 03:22 (unchanged line
numbering checked against a snapshot). All numbers quoted from code outputs were read from
paper/quantum/results/*.json and *.md; independent checks were run in /tmp/claude-0/review4q/.

## Summary

The paper builds a decision-level fault-tolerant cost model for a golden-rule ISC rate written as a QSVT spectral
filter: it checks the factorization of the success probability on a complete CAS(4,4) circuit, costs THC
qubitization of a (12e,11o) Hamiltonian, and compares MLAE, canonical amplitude estimation, direct sampling/SPRT and a
channel-resolved route for two decisions (rate to 20%, triplet-yield class). I re-derived the headline numbers (DKLR,
SPRT, MLAE and Brassard call counts, the walk cost of 3.25e3 Toffolis on 193 qubits with PennyLane 0.45.1, filter
degrees, sampling/MLAE and channel/filter ratios, the quadrature error sum, the lambda sensitivities). They match the
current decision_cost.json and the main-text ranges, so the arithmetic is sound. What I object to is the framing. The
costed problem is solved exactly and classically by the authors' own scripts in minutes. The inputs mix three levels
of theory, and the spin-orbit norm that sets the cost changes 40-fold with the active space. The deepest rigorous
circuits run for years without interruption. The abstract also claims an indifference-zone conclusion that the
procedure does not deliver.

## Recommendation

Major revision. The cost model is carefully done, and most of the fixes are text-only. Two of them (findings 1 and 2)
need new calculations: a scaling or extrapolation to a classically hard, SOC-aware active space, and a cost with a
controlled active space. Without them, the cost numbers cannot support the claim the abstract makes about "what that
decision costs on a fault-tolerant computer".

---

## Major findings

### 1. The costed problem is classically trivial, and the classical baseline is never stated (major)

**Location:** main.tex:40-43 (abstract question), 53-57 (abstract costs), 241-242, 709-713 (Conclusions); main.tex:87-88
cites Babbush et al. for "the classical baseline must be stated".

**What is wrong.** Every Toffoli count is for the (12e,11o) active space. A (12e,11o) CASCI has
C(11,6)^2 = 462^2 = 2.1e5 determinants (M_S = 0), and it is diagonalized exactly on a workstation. The authors do this
themselves: cas_consistent_p.py computes the exact CASCI S1, the T1-T6 energies, the couplings and the golden-rule
k_ISC of the same Hamiltonian. The runs take "seconds": 407 for the BODIPY core, 983 for I2-BDP and 2341 for
aza-BDP-Br2 (results/*.cas_p.json). The decision that the paper costs at 1e11-7e15 Toffolis, years at 1e6 Toffoli/s
(finding 3), is the one the paper already obtains classically in 7-40 minutes. The manuscript never says so. The
introduction cites Babbush et al. on stating the classical baseline, and Section "Fault-tolerant resources, like for
like" compares only with another quantum estimate.

**Evidence.** results/BODIPY_core.cas_p.json, I2-BDP.cas_p.json, azaBDP_Br2.cas_p.json ("exact CASCI", "seconds"
field). main.tex contains no statement that the (12e,11o) problem is classically tractable (grep for "classically",
"tractable", "baseline").

**Fix.**
- (a) Text-only. State in the abstract, in Section "Cost of a screening decision" and in the Conclusions that the
  (12e,11o) numbers are a unit-cost calibration of the method. Say that the same decision is computed exactly and
  classically for this active space in minutes, and that it gives the wrong rate (Section "Active-space
  Hamiltonians").
- (b) New calculation. Give the cost for an active space that is classically hard and contains the orbitals that
  carry the SOC, which is what the Discussion (main.tex:679-681) says is needed. Possible sizes: N = 30-50, or the
  11-45 range of Zhou et al. Give lambda(N) and C_walk(N) from the PennyLane estimator, and alpha_SO for a pi + I 5p /
  C-I sigma* space. Without this, "what that decision costs on a fault-tolerant computer" (abstract) is not
  answered for any problem where a quantum computer would be used.

### 2. The headline costs combine three inconsistent models, and the parameter that sets them varies 40-fold with the active space (major)

**Location:** main.tex:241-242, 572-591, 634-637 (Table 3 caption), 683-689 (Discussion), 53-57 (abstract).

**What is wrong.** The success probability is p = gamma^2 k_ISC/K', with K' proportional to alpha_SO^2. It combines
three inputs from different models:
- k_ISC from the TDDFT/DeltaSCF pipeline;
- alpha_SO from the canonical (12e,11o) window with B3LYP orbitals, a window that the paper shows has 2.2-6.7 of 11
  pi weight and holds iodine orbitals "by accident" (main.tex:590);
- lambda = 17.3 Ha of the BODIPY core, used for every dye.

The cost scales as alpha_SO/V_eff (main.tex:684). For I2-BDP, alpha_SO is 12006 cm^-1 in the costed window
(alpha_so.json) and 299 cm^-1 in the AVAS pi space (Table 2). That is a factor of 40 in the number of amplitude
estimation calls and 1600 in p. The paper argues that neither space is the right one. The abstract and the
Conclusions nevertheless quote the canonical-window numbers ("1e3-1e7 calls", "1e11-7e15 Toffolis") without that
qualification.

**Evidence.** alpha_so.json: I2-BDP 12006, I-BDP 6553, aza-BDP-Br2 855 cm^-1. Table 2 (main.tex:432-444) gives
AVAS 299, 80 and 285 cm^-1 for I2-BDP, I-BDP and I2-BDP-0. With the CASCI energies and couplings of the same window,
p moves by a factor of 1200 down for I2-BDP and 27 up for aza-BDP-Br2 (cas_p.json; main.tex:578-583).

**Fix.**
- (a) Text-only. In the abstract, say that the costs are for a hypothetical (12e,11o) Hamiltonian that reproduces the
  pipeline rates, with the spin-orbit norm of a canonical window, and that the cost scales linearly with
  alpha_SO/V_eff. Present Table 3 as illustrative points on a curve of cost against alpha_SO/V_eff, which
  Figure 2 nearly is already.
- (b) New calculation. Repeat the costs with alpha_SO of a controlled, SOC-aware space (AVAS pi plus iodine 5p or
  C-X sigma*). Report the spread, and use the per-dye lambda where it is available (12.0 Ha for I2-BDP, 11.5 Ha for
  aza-BDP-Br2, thc_regularized.json). This was point Q10 of the previous round and is still open.

### 3. Coherent depth and runtime make the rigorous estimates infeasible, and this is understated (major)

**Location:** main.tex:612-614, 233-237, 598-601; si.tex:186-193.

**What is wrong.** The paper reports total Toffolis and one wall-clock figure ("about 6 years" for D1 of I2-BDP with
MLAE). It does not report the length of a single coherent circuit, which is what fault tolerance has to sustain.
- **Canonical QAE.** The deepest canonical circuit is (2(M-1)+1) x C_call. That is 1.0e14 Toffolis for I2-BDP and
  4.5e14 for I2-BDP-0 (decision_cost.md, "deepest circuit"), so 3.2 and 14 years of uninterrupted coherent
  evolution at 1e6 Toffoli/s. "Needs logical error rates below about 1e-15" (main.tex:614) does not convey that.
  This rigorous estimator is not a usable comparator.
- **MLAE.** The deepest MLAE circuit is not reported. From the schedule in ae_demo.py (100 shots per depth,
  m = 2^k), D1 for I2-BDP needs 9.1e5 calls, so m_max is about 2.3e3. The longest circuit is then about 4.5e3 calls,
  or about 9e11 Toffolis, roughly 10 days. That is a factor of about 100 shorter than canonical. The point favours
  the authors and should be made explicitly.
- **Campaign.** For the 15 dyes the D2 totals alone add up to about 2e15 Toffolis, about 60 years at 1e6 Toffoli/s
  without parallelism.

**Fix.**
- (a) Text and small calculation. Add the MLAE maximum depth to the D1 and D2 tables, and give wall-clock time per
  decision and per campaign. Say that MLAE shots at different depths are independent and can run in parallel on
  several machines.
- (b) New calculation, recommended. Replace canonical QAE as the "rigorous" comparator with a rigorous low-depth
  estimator: iterative QAE (Grinko et al., npj Quantum Inf. 2021), Aaronson-Rall (SOSA 2020) or QSP-based amplitude
  estimation (Rall and Fuller, Quantum 2023). These come with proven confidence intervals, constants close to MLAE
  and no QPE register. Add a physical-qubit and runtime estimate under a stated surface-code model. JCTC readers will
  ask for it.

### 4. "Find that a dye lies within a factor of 2 of the threshold" is not what the D2 procedure delivers (major, because it is in the abstract)

**Location:** main.tex:55-56 (abstract), 610-612, 638-639 (Table 3 footnote a); si.tex:186-193;
decision_cost.py:300-318.

**What is wrong.** The adaptive test stops when the confidence interval excludes one edge of the zone. Its final
angle half-width is dth = max(theta - theta0, theta1 - theta)/2. With an estimate anywhere within ±dth of the truth,
this choice guarantees that the far edge is excluded. It does not guarantee that the near edge is excluded, which
needs dth < (near distance)/2. So the test returns a one-sided answer, for example "k_ISC < 2 k_th", which is a
valid class decision inside the indifference zone. It is not a certificate that the dye lies within a factor of 2 of
the threshold.

**Evidence.** Br2-BDP: p = 1.9e-7, p_th = 2.1e-7 (decision_cost.md, D2 table). This gives theta = 4.36e-4,
theta0 = 3.24e-4 and theta1 = 6.48e-4. Then dth = max(1.12e-4, 2.12e-4)/2 = 1.06e-4, and excluding the near edge
p0 with certainty needs dth < 0.56e-4. For aza-BDP-IPh2, the near distance is 2.8e-5 against 2·dth = 3.8e-5, so the
near edge is not guaranteed either. For a dye close to an edge, a two-sided certificate needs a precision that
diverges.

**Fix.** Text-only: change the abstract and footnote a to "the test returns a class that is correct up to the factor-2
indifference zone". If the authors want the "within a factor of 2" statement, they need a new calculation: cost a
two-sided stopping rule with dth < min(theta - theta0, theta1 - theta)/2 and report it for the two dyes.

### 5. The D1 estimator is incompletely specified: gamma^2 must be measured, and the polynomial bias is not in the 20% (major for correctness, small for cost)

**Location:** main.tex:229-231 ("the estimators see gamma^2 p"), 220-222 (D1 definition); si.tex:172-178, 169;
decision_cost.py:276-279 and docstring "eps = r sqrt(p)/4 so that the induced error in p stays below r p/2".

**What is wrong.**
- **The rate needs gamma^2.** The S1-filter flag is folded into the good state, so amplitude estimation returns
  gamma^2 p. The rate is k_ISC = K' (gamma^2 p)/gamma^2, which needs gamma^2. The paper assumes gamma^2 = 0.5 but
  never measures it, and the model does not cost its estimation. Estimating gamma^2 = 0.5 to about r/2 with MLAE
  takes about 45/(0.1·0.7), roughly 650 calls of C_S1 = 2.6e7, about 2e10 Toffolis. That is negligible against the
  totals, but the estimator is incomplete without it, and its error enters k.
- **The polynomial adds a bias.** The filter precision eps = r sqrt(p)/4 allows a deterministic bias of up to
  r/2 = 10% in p. That comes on top of the 20% statistical error at 95%, so the worst-case D1 error is 30%, not 20%.

**Fix.**
- Text and small calculation: add the gamma^2 estimation (calls and Toffolis) to C_QPE-like fixed costs.
- Either state D1 as "20% statistical plus at most 10% polynomial bias" or split the 20% between the two.
- The CAS(4,4) check uses the exact gamma^2 = 0.94 (full_circuit_cas44.json). Say so.

---

## Minor findings

### 6. The compiled si.pdf is stale against si_tables.tex, and the S9 caption disagrees with Table 3 (minor, text-only and recompile)

**Location:** si.pdf Table S9; si_tables.tex:70-106; main.tex:636-639.

**Evidence.** si.pdf (built 01:10) still shows the oracle-margin D2 of an earlier model:

| Dye | D2, si.pdf | D2, main Table 3 and si_tables.tex |
|---|---|---|
| BDP | 9.6e10 | 3.1e11 |
| I2-BDP | 1.8e13 | 8.8e13 |
| Br2-BDP | 2.8e14 | 9.5e13 |
| aza-BDP-Br2 | 6.4e11 | 3.3e12 |

si.pdf also has older D1 values (I2-BDP MLAE 1.6e14 against 1.8e14). BDP at 9.6e10 contradicts the main text's
"1e11-2e12" for light-atom negatives. si_tables.tex (regenerated 03:22) now matches decision_cost.json.

**Fix.** Recompile si.pdf. The S9 caption says "D2 ... by MLAE", while Table 3 and the text use "best estimator".
For aza-BDP-Ph4 and BDP-0 the best estimator is sampling (1.3e11 and 1.6e11 against MLAE 1.6e11 and 1.8e11). Make
the two consistent, or add a sampling column.

### 7. The definition of p changes between eq 3, the text and Figure 2 (minor, text and figure label)

**Location:** main.tex:215 (eq 3, p = k/K'); main.tex:576, 636, 664 (p = gamma^2 k/K'); decision_cost.py:489 (axis
label "p = k_ISC/K'" while the plotted points are r["p"], which includes gamma^2).

**Fix.** Use one symbol, for example p and p_good = gamma^2 p, and relabel the axis.

### 8. The bare-Hamiltonian comparison omits the core (minor, text-only)

**Location:** main.tex:583, "1200 times smaller or 27 times larger".

**Evidence.** The BODIPY core has p_CAS = 6.4e-13 (BODIPY_core.cas_p.json) against 7.8e-16 in the pipeline, which
is 815 times larger.

**Fix.** Say "800 times larger for the core, 1200 times smaller for I2-BDP and 27 times larger for aza-BDP-Br2".

### 9. Walk steps per call are understated (minor, text-only)

**Location:** main.tex:553-555, "about 3e4 walk steps ... about 1e8 Toffolis".

**Evidence.** decision_cost.json gives degrees of 3.8-5.9e4 for D1 and 2.1-5.5e4 for D2, so the filter costs
0.7-1.9e8 Toffolis.

**Fix.** Say "2-6e4 walk steps (0.7-1.9e8 Toffolis)". This was minor point 3 of the previous round and is still open.

### 10. The split of 1/sqrt(p) into alpha_SO/||chi|| and ||chi||/V_eff mixes levels of theory (minor)

**Location:** main.tex:584-586, 685-689.

**What is wrong.** ||chi|| = 1120 cm^-1 is a TDDFT quantity (decision_cost.md), and it is essentially the single S1-T5
coupling of 1119 cm^-1 at 3.80 eV (Table S6). The CASCI chi norm of the costed Hamiltonian over T1-T6 is 114 cm^-1
(I2-BDP.cas_p.json). In a consistent model alpha_SO/||chi|| is at most about 105, not 11.

**Fix.** Small calculation: compute ||chi||^2 = <S1|H_SO^2|S1> exactly in the CASCI (it is a one-step
contraction). Alternatively, drop the split and keep alpha_SO/V_eff. Previous Q3 still partly open.

### 11. The S1-preparation assumptions are not supported by the costed Hamiltonian (minor, small calculation)

**Location:** si.tex:172-178, decision_cost.py:84-87.

**What is wrong.**
- The S1 filter assumes a 0.5 eV S1-S2 gap. The CASCI of the costed window gives 0.34 eV for aza-BDP-Br2, the top
  candidate (E_S2 - E_S1 = 3.878 - 3.534 eV, azaBDP_Br2.cas_p.json). That raises C_S1 by 1.5 times.
- gamma^2 = 0.5 and D = 1e4 CSFs are assumptions. gamma^2 could be taken from the CASCI S1 vector (weight of the
  leading 1e4 CSFs).
- The sum-of-Slaters formula c_prep is for determinants. A CSF expansion needs a few determinants per CSF, which is
  negligible but should be stated.

**Fix.**
- Take the gap and gamma^2 from the CASCI.
- State the leakage bias: S0 or S2 components of amplitude eps_S = 1e-3 pass through H_SO and the filter with
  weight about eps_S^2 (V_other/V_eff)^2. That is negligible for the heavy-atom D1, but not small for the
  light-atom dyes, where V_eff is 1e-5-1e-2 cm^-1.

### 12. Qubit and walk-cost bookkeeping (minor, text and rerun)

**Location:** main.tex:555-558; logical_qubits.py; thc_pennylane_resources.md; decision_cost.md header.

**What is wrong.**
- logical_qubits.py uses lambda = 15.9 Ha rather than the 17.3 Ha of the cost model, and a stale filter degree of
  51264 (the current D1 maximum is 58685). I checked that the wire count is unaffected: 193 at both lambdas, rank 33,
  5 batches.
- It computes 19 phase bits. The canonical model needs M = 2^20 for I2-BDP-0, so 20 bits and 221 qubits in total.
- thc_pennylane_resources.md lists ranks 22 and 44 only, while C_WALK = 3.25e3 is hard-coded. I reproduced
  3.251e3 Toffolis and 193 wires at rank 33 with PennyLane 0.45.1; add that row.
- decision_cost.md still says the degree is "fitted", while the code interpolates the direct grid.
- Registers for the trial-state QROM and the H_SO coefficients are excluded from the qubit count, although their
  Toffolis are counted.

**Fix.** Rerun the scripts and state the exclusions with an estimate.

### 13. C_SO counts spatial rather than spin orbitals (minor, text-only)

**Location:** decision_cost.py c_so(n=11).

**What is wrong.** The SOC matrix is a complex 22x22 spin-orbital matrix. Diagonalizing it takes 231 complex Givens
rotations, applied and undone, so about 2-4e4 Toffolis rather than 4.5e3. That is still negligible. The
logical_qubits.py docstring already says "22 spin orbitals". Previous minor 11 is still open.

### 14. The MLAE constant and its range of use (minor, text-only)

**Location:** si.tex:182-185, 200-204; ae_demo.json.

**What is wrong.**
- c = 45 is the median of c_eff = 35.7, 43.7, 45.0, 39.8, 56.2, 106 and 255 over K = 7-13. That is not a constant:
  for K >= 11 the implementation stalls (q95 at the true theta is 0.087-0.094).
- I checked the Cramér-Rao bound at K = 8 (Fisher information 4·100·Σ(2m+1)^2). It gives a relative sigma of
  0.0715 and a two-sided q95 of 0.140, against the simulated 0.146. So the constant is right near the 20% target.
  The stall is an artefact of the hierarchical grid (w = 30/sqrt(F)), not of MLAE.
- For D2 on the light-atom negatives the constant is applied at theta close to 0 and relative errors of 1e5 or
  more, a regime that was not simulated.

**Fix.** Quote the Cramér-Rao-based constant, which is θ-independent because the Fisher information is, as the
justification. Say explicitly where the calibration was extrapolated.

### 15. The D2 sampling and MLAE comparison is not like for like (minor, text-only)

**Location:** si.tex:182-193; main.tex:604-608.

**What is wrong.**
- SPRT is costed by its expected sample number from Wald's approximations, while MLAE is costed by a fixed
  95th-percentile budget. This favours sampling, so the "4-270 times" claim survives, but say so.
- The per-stage confidence 0.05/L uses L computed from the true theta. A practical test needs a schedule fixed in
  advance, for example delta_k = 6 delta/(pi^2 k^2). The effect on cost is small.

### 16. Error terms missing from the near-resonance budget (minor, text-only)

**Location:** main.tex:566-570.

**What is wrong.**
- The quadrature sum (0.003, 0.027 and 0.017 eV, giving 0.032 eV; I verified it) omits the block-encoding precision.
  The coefficient and rotation bits (20 and 27) are set for 0.16 mHa, about 0.004 eV.
- The degree depends on the single-mode MLJ width: lambda_s = 0.05 eV raises the degree, and the D1 cost, by
  sqrt(2). A smoother multimode line shape (SI Section S12) would lower it.
- The QPE cost pi lambda'/(2 eps) is an RMS-type precision, and 95% confidence needs a small constant factor
  (previous minor 13).

**Fix.** State all three.

### 17. Channel-resolved route: favourable assumptions not all stated (minor, text-only)

**Location:** main.tex:616-624; decision_cost.py:234-255.

**What is wrong.**
- The code charges each T_n phase estimation as c_qpe_s1. That implicitly assumes a separate T_n trial state with
  gamma^2 = 0.5.
- The T_n filter width is the distance to the nearest triplet in the TDDFT list. The CAS triplet spectrum is denser:
  in I2-BDP CASCI, T3, T4 and T5 lie within 0.033 eV of each other (I2-BDP.cas_p.json). That does not matter for the
  T1 channels costed here, but it does for upper-triplet channels.

**Fix.** Add both points to the existing caveat.

### 18. The abstract gives optimistic numbers without saying so (minor, text-only)

**Location:** main.tex:53-57.

**What is wrong.**
- "Calls of about 2e8 Toffolis" should be 0.9-2.2e8.
- The quoted totals are MLAE or best-estimator values. MLAE is called "the optimistic estimate" at main.tex:235.
  The rigorous canonical totals, 5-9 times higher for D1, do not appear in the abstract.

**Fix.** Add "with maximum-likelihood amplitude estimation; rigorous bounds are 5-10 times higher", and the
qualification of finding 2.

### 19. Missing prior work (minor, text-only)

**Location:** Introduction (main.tex:84-88), Methods (main.tex:226-238), Discussion.

Please cite and, where relevant, compare with the following:
- **Eigenstate filtering** for the S1 preparation: Lin and Tong, Quantum 4, 372 (2020).
- **QSVT and phase factors** at degree 5e4: Martyn et al., PRX Quantum 2021; Dong, Meng, Whaley and Lin, PRA 2021.
- **Rigorous low-depth amplitude estimation**: Grinko et al. 2021; Aaronson and Rall 2020; Rall and Fuller 2023;
  Giurgica-Tiron et al., PRR 2022.
- **Heisenberg-limited expectation-value estimation**: Knill, Ortiz and Somma, PRA 2007; Huggins et al., PRL 2022.
- **A fault-tolerant spectral observable through an operator block encoding**: Fomichev et al., X-ray absorption
  spectroscopy, 2024. Its dipole-normalization issue is the analogue of alpha_SO/||chi||.
- **Campaign- and decision-level feasibility studies**: Bellonzi et al. on homogeneous catalyst discovery, 2024;
  Goings et al., PNAS 2022.
- **One-norm reduction**, which acts on the same "cost lever": BLISS (Loaiza and Izmaylov, JCTC 2023) and
  symmetry-compressed double factorization (Rocca et al., JCTC 2024).

---

## Checks that passed (for the editor)

- **D1 call counts.** DKLR Y1 = 319, giving 5.2e9 calls for I2-BDP. MLAE 45/(0.2 sqrt(6.1e-8)) = 9.1e5. Brassard
  M = 2^18 with 16 runs gives 8.4e6. All match decision_cost.md.
- **Main-text ratios.** Sampling/MLAE for D1 is 301-25423, quoted as "3e2-3e4". Channel/filter is 11.6-27.8,
  quoted as "12-28". D2 sampling/MLAE for heavy-atom dyes is 3.6-274, quoted as "4-270". Canonical/sampling is
  1.24-9.4, quoted as "1.2-10". The cost of not knowing the margin is 2.2-5.3, quoted as "2-5".
- **D2 SPRT, I2-BDP.** By hand, ln(19)/(p ln 4 - p1 + p0) = 3.7e7, as reported.
- **Lambda sensitivity.** lambda = 15.9 lowers the costs by 7% and 37.5 doubles them, consistent with the
  _lambda*.md files.
- **Walk cost.** PennyLane 0.45.1, rank 33, 5 batches: 3.251e3 Toffolis and 193 wires, and the walk is Clifford plus
  Toffoli only.
- **CAS(4,4) full circuit.** p_call/p_pred = 1.0041, and gamma^2 = 0.94.
