# Related work and gap analysis

Literature review for the psscreen paper (JCTC target). Search done on 2026-11-16. Every work cited here
was checked against its Crossref record (DOI) or its arXiv listing; BibTeX entries are in
`paper/related_work.bib` (75 entries, keys do not collide with `paper/refs.bib`). Works that could not be
verified are listed at the end and are **not** in the .bib file.

Our four claims, as used throughout this document:

- **C1. Precision budget.** S–T gaps must be known to about 0.07 eV for k_ISC within a factor of 2 (MLJ
  model), with absorption-window edge sensitivity and a fragile top candidate (tolerance 0.016 eV).
- **C2. Gap-resolved ISC observable.** Fermi's golden rule written as a QSVT spectral filter on
  H_SO|S1⟩, verified in PennyLane, at about the cost of the arXiv:2512.15889 proxy (~10⁸ Toffolis per shot,
  ~130–190 logical qubits with THC).
- **C3. Dynamic correlation is the bottleneck.** With CASCI energies the correct observable mis-ranks the
  dyes (S1 off by 1.4–2.0 eV at CAS(12,11)); the useful quantum target is energies with dynamic correlation.
- **C4. Experimental benchmark** on nine BODIPY dyes with measured absorption, k_ISC and Φ_T.

---

## 1. Bottom line

1. **No work found sets an end-to-end accuracy requirement for a quantum-computed photophysical
   observable from the screening decision it serves (C1).** Quantum resource papers fix the precision a
   priori: chemical accuracy (1.6 mHa) for energies, or an additive sampling error (10% in
   arXiv:2512.15889), or a broadening matched to experimental resolution (≈1 eV for XAS,
   Fomichev2024XAS). Classical ISC papers note that rates depend steeply on ΔE_ST (Marian2012,
   Penfold2018, Marian2021) and Shi2022 tests how parameter accuracy propagates into high-throughput
   non-radiative rates, but none of them turns that into a tolerance for a quantum algorithm. C1 is the
   cleanest novelty claim. It should be phrased as "decision-linked tolerance", not as "the first
   sensitivity analysis of ISC rates", because classical sensitivity analyses exist.
2. **C2 is a new application, not a new algorithm.** Estimating ⟨ψ|f(H)|ψ⟩ for a smooth f with
   QSP/QSVT or an integral transform is established: Roggero2019 (linear response via integral
   transforms), Roggero2020 (Gaussian integral transform), Rall2020 (spectral densities and related
   quantities with block encodings), Lu2021 (energy filters), Kosugi2020 (optical absorption), and, inside
   chemistry, Fomichev2024XAS and Loaiza2025 (spectroscopy via QPE). Because the MLJ line shape is a
   Poisson-weighted sum of Gaussians, our filter is literally a sum of Gaussian integral transforms in the
   sense of Roggero2020. The golden-rule ISC rate has been *written down* for a fault-tolerant quantum
   workflow before (Zak2026, as a simplified k ∝ |V|²/ΔE² and as a short-time proxy, with no
   implementation or numbers), and arXiv:2512.15889 states the golden rule as background but implements
   only the gap-free proxy. What we found nobody has done: implement the golden-rule rate over all triplets
   as one QSVT filter on H_SO|S1⟩, verify it numerically, and cost it against the proxy. Claim exactly that
   and cite the spectral-density literature as the method's origin.
3. **C3 is well known in general and new in its specific consequence.** That a bare active space misses
   BODIPY excitations by ~1–2 eV and needs CASPT2/NEVPT2 is established classically (Momeni2015,
   Chibani2013, BerraudPache2020). That quantum active-space solvers need a dynamic-correlation layer is a
   recognized research line (Takeshita2020, Krompiec2022, Huang2023DSRG, Motta2020TC, Kumar2022TC,
   Mitarai2022). What is new is the demonstration that, for the PDT observables proposed for fault-tolerant
   hardware, an *exact* active-space solution ranks dyes in the wrong order, so the quantum resource
   estimates in the literature cost the wrong problem. Present C3 as a quantitative consequence of a known
   limitation, with the literature above cited up front.
4. **C4 has partial precedents.** Nykanen2024 benchmarks ΔUCCSD/ΔADAPT-VQE S1 energies of six BODIPY
   dyes against experiment (MAE 0.13 eV for ΔUCCSD(8,8) in simulation), and Genin2025 benchmarks
   classically emulated iQCC T1 energies of Ir/Pt emitters against experiment (MAE 0.05 eV). Neither
   benchmarks ISC rates or triplet yields, and neither compares a fault-tolerant observable with a
   classical pipeline on the same dyes. Our benchmark is the first we found that tests a proposed quantum
   ISC observable against measured k_ISC and Φ_T.

### Novelty risks (read before finalizing the introduction)

| Risk | Work | What overlaps | How to handle it |
|---|---|---|---|
| **Medium** | Zak2026 (arXiv:2601.18898, Jan 2026) | Fault-tolerant workflow for PDT photosensitizers (BODIPY named), including SOC matrix elements and a golden-rule ISC rate with a vibronic density of states, plus a short-time population proxy. | Cite explicitly. Differences: Zak gives no implementation, no numbers and no gap tolerance; the rate is reduced to k ∝ \|⟨S1\|H_SO\|T1⟩\|²/ΔE² for S1→T1 only; the algorithm is QPE on the Pauli-Breit Hamiltonian, not a spectral filter on H_SO\|S1⟩. |
| **Medium** | Rall2020, Roggero2019/2020, Lu2021 | The mathematical primitive behind C2 (spectral density / integral transform / energy filter via block encodings). | Say "we apply the spectral-density estimation framework of refs X to the golden-rule ISC rate". Do not claim a new algorithmic primitive. Note the Gaussian-transform link to the MLJ line shape. |
| **Low–medium** | Nykanen2024 (Algorithmiq) | Quantum-algorithm excitation energies of BODIPY compared with experiment, motivated by PDT. | Cite in the C4 paragraph. Ours differs in property (ISC, Φ_T), in method class (fault-tolerant observables vs VQE) and in the head-to-head with a classical pipeline. |
| **Low** | Genin2025 (OTI Lumionics / Samsung) | "Fault-tolerant-designed" solver benchmarked against experiment for triplet energies; argues classical tractability up to ~200 qubits. | Cite as a parallel benchmark-against-experiment effort in OLED phosphors. |
| **Low** | Algorithmiq Q4Bio perspective (white paper, no DOI or arXiv id) | Frames PDT as a driver for quantum advantage; names ISC pathways and singlet–triplet energetics as goals. | Grey literature; mention only if a reviewer raises it. Not in the .bib. |
| **Low** | Baig2024 | Ab initio ISC dynamics of iodo-BODIPY with a golden-rule α⁻² extrapolation of SOC. | Classical; cite as a contrasting high-cost route to k_ISC (their k_ISC ≈ 1.7×10¹¹ s⁻¹, Φ_T 0.85). |

**Housekeeping:** arXiv:2512.15889 is now published in Phys. Rev. A (DOI 10.1103/nj6b-msgw, issued
2026-11-04; entry `Zhou2026PRA`). The paper should cite the journal version and re-check our description
of it against that version, since the published text may differ from v1.

---

## 2. Gap table

Legend: ✓ does it; ◐ partially or qualitatively; ✗ does not. "Exp." = compared with experiment.

| Work | Domain | C1 decision-linked accuracy | C2 gap-resolved ISC on QC | C3 dynamic-corr. bottleneck for the observable | C4 exp. benchmark | Note |
|---|---|---|---|---|---|---|
| Zhou2026PRA / arXiv:2512.15889 | FT QC, PDT | ✗ (fixed ε_samp = 10%) | ✗ (gap-free S1→T1 proxy; golden rule only as background) | ✗ (CASPT2 cited for the active space, not in the workflow) | ✗ | The work we respond to |
| Zak2026 | FT QC, PDT | ✗ | ◐ (golden rule written, k∝V²/ΔE², not implemented) | ✗ | ✗ | Closest overlap with C2 |
| Nykanen2024 | VQE, BODIPY | ✗ | ✗ | ◐ (ΔUCCSD error falls with active space) | ✓ S1 only, 6 dyes | Closest overlap with C4 |
| Genin2025 | iQCC, OLED | ✗ | ✗ | ✗ | ✓ T1 only | Emulated FT solver |
| Gao2021 | VQE-type, TADF | ✗ | ✗ (S1, T1 energies only) | ✗ | not checked | Hardware demo |
| Fomichev2024XAS | FT QC, spectroscopy | ◐ (1 eV broadening = exp. resolution) | ✗ (XAS, not ISC) | ✗ | ✗ | Same primitive family |
| Loaiza2025 | FT QC, nonlinear spectroscopy | ✗ | ✗ | ✗ | ✗ | Generalized QPE |
| Roggero2019, Roggero2020, Rall2020, Lu2021 | QC algorithms | ✗ | ◐ (primitive) | ✗ | ✗ | Origin of C2's method |
| Kosugi2020 | QC linear response | ✗ | ✗ | ✗ | ✗ | Optical absorption |
| Sawaya2019, Huh2015 | QC vibronic spectra | ✗ | ◐ (FC factors, not ISC) | ✗ | ✗ | Line-shape side of C2 |
| Motlagh2024, Pocrnic2026, Ollitrault2020NAD | QC dynamics | ✗ | ◐ (could give ISC from dynamics) | ✗ | ✗ | Higher cost, beyond MLJ |
| Reiher2017, vonBurg2021, Lee2021, Goings2022, Kim2022, Rocca2024, Caesura2025 | Resource estimation | ✗ (chemical accuracy assumed) | ✗ | ◐ (Goings2022 checks classical solvability) | ✗ | Cost models we reuse |
| Steudtner2023 | FT observables | ◐ (precision for observables) | ✗ | ✗ | ✗ | Cost of expectation values |
| LeeS2023, Chan2024, Liu2022, Elfving2020, Babbush2021 | Utility debates | ◐ (general) | ✗ | ◐ (classical heuristics argument) | ✗ | Framing |
| Takeshita2020, Krompiec2022, Huang2023DSRG, Motta2020TC, Kumar2022TC, Mitarai2022 | Dyn. corr. on QC | ✗ | ✗ | ◐ (recognize the need; ground states/small molecules) | ✗ | Directions for C3 |
| Momeni2015, Chibani2013, BerraudPache2020, Jacquemin2012 | Classical BODIPY/cyanine | ✗ | n/a | ✓ (classical: CAS needs PT2; TDDFT cyanine error) | ✓ absorption | Classical basis for C3 |
| Marian2012/2021, Penfold2018, Etinski2011, deSouza2019 | Classical ISC theory | ◐ (gap dependence qualitative) | n/a | ✗ | ◐ | Golden rule, spin-vibronic |
| Shi2022 | Classical HT rates | ◐ (parameter accuracy vs rate) | n/a | ✗ | ◐ trends | Closest classical analogue of C1 |
| Baig2024 | Classical ISC dynamics, I-BODIPY | ✗ | n/a | ✗ | ◐ Φ_T only | k_ISC 1.7×10¹¹ s⁻¹ |
| Ly2021 | Experiment | n/a | n/a | n/a | data | Our reference data |
| Loos2023 | Classical ST gaps | ◐ (chemically accurate ΔE_ST) | n/a | ✓ (needs high-level CC) | ✗ | ST-gap accuracy context |
| GomezBombarelli2016, Chen2025, Fallani2026, Westermayr2021 | HT / ML screening | ✗ | n/a | ✗ | ◐ (GB2016 synthesized) | Screening context |
| **This work** | | ✓ | ✓ | ✓ | ✓ (9 dyes, λ_max, k_ISC, Φ_T) | |

---

## 3. Review by area

### 3.1 Quantum computing for photosensitizers and PDT

- **Zhou et al.** (`Zhou2026PRA`, arXiv:2512.15889). Fault-tolerant algorithms for PDT screening: QSVT
  threshold projection of the dipole spectral weight in a therapeutic window and a short-time
  |⟨T1,M|e^{-iH_SO t}|S1⟩|² ISC proxy, with BLISS-THC resource estimates (Orig-BDP CAS(12,11): 177 qubits,
  2.7×10⁷ Toffolis for the window, 7.96×10⁷ for the proxy; up to ~350 qubits and 2.5–5.8×10⁹ Toffolis at
  45 orbitals). The text states that the proxy "is not an absolute ISC rate", that k_ISC "is governed by
  Fermi's golden rule", and that a vibronic route is needed in the weak-SOC regime. No ISC values, no
  spectra and no comparison with experiment are reported (details in `xanadu_comparison.md`).
- **Zak** (`Zak2026`, arXiv:2601.18898). Block encoding of the double-factorized Pauli–Breit Hamiltonian
  with spin-controlled SWAP networks (~2× lower LCU prefactor), embedded in a hybrid PDT design workflow
  (BODIPY, aza-BODIPY, phthalocyanines, Pt-BODIPY). SOC matrix elements by Hadamard test or amplitude
  estimation; ISC as a golden-rule expression with a vibronic density of states, used in simplified form
  k_ISC ∝ |⟨S1|H_SO|T1⟩|²/ΔE², and a short-time population proxy. Only asymptotic costs, no computed
  molecules, no experiment, no gap-accuracy discussion. **Must be cited**; see the novelty table.
- **Nykänen et al.** (`Nykanen2024`, arXiv:2404.16149, Algorithmiq / Cleveland Clinic). ΔUCCSD-VQE and
  ΔADAPT-VQE (non-Aufbau ΔSCF references) for S1 of six BODIPY derivatives, statevector simulation with
  CPCM. MAE vs experiment: ΔCCSD 0.054 eV, ΔUCCSD(8,8) 0.131 eV, TDDFT/B3LYP 0.368 eV, EOM-CCSD 0.498 eV.
  No triplets, SOC or ISC. Useful contrast: our ΔSCF-corrected B3LYP pipeline reaches a mean error of
  −0.055 eV on compact dyes classically, at the same level as their ΔCCSD.
- **Algorithmiq Q4Bio perspective** ("Photodynamic Therapy as a Driver for Useful Quantum Advantage in
  Chemistry", white paper on algorithmiq.fi, no DOI or arXiv id; not in the .bib). Near-term pipeline
  (ΔADAPT-VQE, informationally complete measurements, quantum-boosted DMRG) for TLD1433 on IBM hardware;
  ISC pathways and S–T energetics are named as future goals; no fault-tolerant estimate. Grey literature.
- **Genin et al.** (`Genin2025`, arXiv:2512.13657v2, OTI Lumionics / Samsung). iQCC, a fault-tolerant-
  oriented solver, run on classical hardware at ~200 logical qubits for T1 energies of Ir(III)/Pt(II)
  phosphors; MAE 0.05 eV and R² 0.94 vs 77 K emission. Argues these systems stay classically tractable.
- **Gao et al.** (`Gao2021`, IBM/Mitsubishi Chemical/Keio). Variational quantum calculations of S1 and T1
  of phenylsulfonyl-carbazole TADF emitters, including runs on IBM hardware; the target quantity is ΔE_ST, the same
  quantity whose required accuracy we quantify.

### 3.2 Quantum algorithms for excited states, spectra and photochemical dynamics

- Near-term excited-state methods: QSE (`McClean2017`), VQD (`Higgott2019`), SSVQE (`Nakanishi2019`),
  MC-VQE (`Parrish2019`), qEOM (`Ollitrault2020qEOM`). All deliver state energies; none delivers a rate.
- Fault-tolerant spectroscopy: XAS via frequency-domain, QPE-sampling and time-domain methods, with core-
  valence separation or filtering and a 1 eV target broadening (`Fomichev2024XAS`, ~900 logical qubits and
  ~10⁸ Toffolis per run for CAS(22e,18o)); nonlinear spectroscopy via generalized QPE (`Loaiza2025`).
- Spectral densities, linear response and filters: `Roggero2019` (integral-transform linear response),
  `Roggero2020` (Gaussian integral transform), `Rall2020` (estimating ⟨ψ|f(H)|ψ⟩, densities of states and
  related quantities from block encodings), `Lu2021` (energy filters for finite-energy states),
  `Kosugi2020` (molecular optical absorption), `Sakuma2025` (QPE-based filters, Kaiser windows, DOS).
  **These are the methodological ancestors of C2.**
- Vibronic and nonadiabatic: `Huh2015` (boson sampling for Franck–Condon profiles), `Sawaya2019` (qubit
  algorithm for vibronic spectra), `Ollitrault2020NAD` (nonadiabatic dynamics), `Motlagh2024` (vibronic
  Hamiltonians with many states and modes, singlet fission), `Pocrnic2026` (pre-Born–Oppenheimer
  dynamics, ~1.4×10³ qubits for NH₃+BF₃). These could supply the Franck–Condon factor that our MLJ model
  treats classically, at much higher cost.
- Rate constants on a quantum computer: `Lidar1999` (thermal rate constants via flux correlation, an
  early precedent for computing rates rather than energies).
- Primitives: QSVT (`Gilyen2019`, `Martyn2021`), qubitization (`Low2019`), ground-state preparation
  (`LinTong2020`), initial-state preparation for chemistry (`Fomichev2024ISP`), cost of observables
  (`Steudtner2023`).

### 3.3 Fault-tolerant resource estimation

`Reiher2017` (FeMoco, Trotter QPE), `Babbush2018` (linear-T encodings), `Berry2019` (single and double
factorization), `vonBurg2021` (double factorization, catalysis), `Lee2021` (THC), `Goings2022` (P450, with a
careful classical-solvability check), `Kim2022` (Li-ion electrolytes, PsiQuantum), `Loaiza2023BLISS` (BLISS),
`Rocca2024` (symmetry-compressed DF), `Caesura2025` (BLISS-THC with active-volume compilation: P450
(63e,58o) λ 389→131 Ha, 7.8×10⁹→1.7×10⁹ Toffolis). All assume a target precision (usually 1.6 mHa) rather
than derive one from the downstream decision. Our THC numbers reproduce the arXiv:2512.15889 qubit count
within the angle-loading trade-off (`xanadu_comparison.md`), so this literature is where we cite the cost
model, not where we claim novelty.

### 3.4 Accuracy requirements and utility of quantum chemistry on quantum computers

Reviews: `Cao2019`, `McArdle2020`, `Bauer2020`, `Liu2022`. Utility and advantage debates: `Babbush2021`
(focus beyond quadratic speedups), `Elfving2020` (industrial relevance), `LeeS2023` (evidence for
exponential advantage in ground-state chemistry is lacking), `Chan2024` (classical heuristics),
`Santagati2024` and `Blunt2022` (drug discovery). The recurring point is that "chemical accuracy in the
active space" is not the same as accuracy with respect to experiment. None of these sets a property-
specific tolerance for a photophysical screening decision; C1 adds one.

### 3.5 Dynamic correlation on top of quantum active-space solvers

`Takeshita2020` (virtual-space correction from RDMs), `Krompiec2022` (SC-NEVPT2 from quantum-measured RDMs,
4-RDM cumulant issues), `Huang2023DSRG` (DSRG-downfolded Hamiltonians), `Motta2020TC` and `Kumar2022TC`
(transcorrelated Hamiltonians, including excited states in Kumar2022TC), `Mitarai2022` (perturbation
theory with QSP). These target ground states or small molecules; none evaluates whether the corrected
energies reach the ~0.07 eV gap tolerance on dye-sized systems. This is the natural "where quantum could
still help" paragraph for C3.

### 3.6 Classical ISC theory and computation

Reviews: `Marian2012`, `Marian2021`, `Penfold2018` (spin-vibronic mechanism), `Baryshnikov2017`
(phosphorescence). Rate methods: `Etinski2011` (time-dependent golden-rule ISC rates), `deSouza2019`
(path-integral excited-state dynamics in ORCA), `Shi2022` (reduced-mode, high-throughput IC/ISC rates with
a parameter-accuracy analysis), MLJ (`Jortner1976`, in refs.bib). El-Sayed rules (`ElSayed1963`, in
refs.bib). The steep ΔE dependence of k_ISC that C1 quantifies is qualitatively textbook here; the
quantitative tolerance against a screening decision is not.

### 3.7 BODIPY photophysics and the TDDFT/cyanine problem

PDT and BODIPY background (`Agostinis2011`, `Loudet2007`, `Kamkaew2013` in refs.bib; `Zhao2013` for
triplet photosensitizers). Benchmarks: `Chibani2013` (BODIPY optical signatures with ab initio tools),
`LeGuennic2015` (cyanine challenge), `Jacquemin2012` (TDDFT "not guilty" for cyanines, a nuance worth
citing alongside), `Momeni2015` (TDDFT vs CC2/CASPT2/NEVPT2 for BODIPY/aza-BODIPY; also the source of the
CAS(12,11) active space used by arXiv:2512.15889), `Momeni2016` (heavy-atom-free BODIPY dimers),
`BerraudPache2020` (bt-PNO-STEOM-CCSD for BODIPY). Experimental reference: `Ly2021` (in refs.bib). ISC
dynamics: `Baig2024`.

### 3.8 ΔSCF and orbital-optimized excited states

`Ziegler1977`, `Gilbert2008` (in refs.bib), `Kowalczyk2011` (ΔSCF for organic dyes, the closest precedent
for our ΔSCF correction), `Hait2021` (orbital-optimized DFT review). Nykanen2024 also builds on ΔSCF
references, which is a useful bridge sentence.

### 3.9 TADF, singlet–triplet gaps and high-throughput / ML screening

`Uoyama2012` (TADF), `Loos2023` (chemically accurate inverted ΔE_ST; shows how much method quality the
sign and size of a gap demand), `GomezBombarelli2016` (high-throughput virtual screening with experimental
validation), `Westermayr2021` (ML for excited states), `Chen2025` (generative deep learning for heavy-atom-
free triplet photosensitizers using ΔE_ST and absorption filters from TDDFT labels), `Fallani2026` (active
learning for Type I transition-metal PDT photosensitizers). These screens use gap thresholds of a few
tenths of an eV (e.g. ΔE_ST ≤ 0.30 eV in Chen2025) without an analysis of how label errors change the
selected set, which is the classical counterpart of C1.

---

## 4. Suggested text for the paper

A draft "Related work" paragraph the paper thread can adapt (citation keys from both .bib files):

> Quantum algorithms for excited states range from near-term variational methods
> \cite{McClean2017,Higgott2019,Nakanishi2019,Parrish2019,Ollitrault2020qEOM} to fault-tolerant
> spectroscopy based on phase estimation and spectral filters
> \cite{Fomichev2024XAS,Loaiza2025,Kosugi2020}. Estimating a smooth function of the Hamiltonian in a
> given state, the primitive behind spectral densities and linear response, has been formulated with
> integral transforms and block encodings \cite{Roggero2019,Roggero2020,Rall2020,Lu2021}. For
> photodynamic therapy, Zhou \emph{et al.}\ proposed a QSVT window projection and a short-time ISC proxy
> with THC resource estimates \cite{Zhou2026PRA}, Zak outlined a Pauli--Breit workflow in which ISC enters
> through a golden-rule expression \cite{Zak2026}, and Nyk\"anen \emph{et al.}\ benchmarked variational
> excitation energies of BODIPY dyes against experiment \cite{Nykanen2024}. Resource estimates for
> chemistry typically assume chemical accuracy in the active space
> \cite{Reiher2017,vonBurg2021,Lee2021,Goings2022,Caesura2025}, while it is well established that bare
> active spaces miss the excitation energies of BODIPY and cyanine dyes without dynamic correlation
> \cite{Momeni2015,Chibani2013,BerraudPache2020}, and several schemes add dynamic correlation to quantum
> active-space solvers \cite{Takeshita2020,Krompiec2022,Huang2023DSRG,Kumar2022TC}. What has not been
> done is to ask what accuracy the screening decision itself requires, to evaluate a mechanism-resolved
> ISC observable at the cost of the proposed proxy, and to test both against measured ISC rates. We do
> this here, writing the golden-rule rate as a spectral filter in the sense of
> refs~\citenum{Roggero2020,Rall2020}.

---

## 5. Verification log

Verified against Crossref (DOI record): all `@article`/`@inproceedings` entries in `related_work.bib`.
Verified against arXiv/alphaXiv listing (title, authors, id): Zak2026, Nykanen2024, Genin2025,
Motlagh2024, Pocrnic2026, Fomichev2024XAS, Sakuma2025, Caesura2025, Elfving2020, Krompiec2022,
Mitarai2022, Baig2024, Fallani2026. Open items marked in the .bib `note` fields: Zhou2026PRA volume and
article number (not yet in Crossref), Mitarai2022 journal version, deSouza2019 volume/pages.

Considered and **left out** (could not be verified in this pass, or off-topic): Send/Valsson/Filippi
cyanine QMC benchmark (JCTC 2011), Pollice et al. inverted-gap molecules (Matter 2021), Samanta et al.
RISC rates (JACS 2017), Lower & El-Sayed review (Chem. Rev. 1966), Bauman et al. downfolding (JCP 2019),
Kassal et al. chemical dynamics (PNAS 2008), Tchapet Njafa et al. TADF active learning (arXiv:2512.06029,
TADF-only), Rossi two-electron PDT model (arXiv:2301.03653, pedagogical). Add any of these only after
checking the DOI.
