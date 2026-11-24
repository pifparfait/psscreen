# Pre-review read (associate editor) and third report (broad theoretical chemist)

Manuscript: "How Accurate Must Excited States Be to Screen Triplet Photosensitizers?" (main.tex / main.pdf, 32 pp.; si.tex / si.pdf, 21 pp.). Checked against the compiled PDFs (pdftotext), main.log/si.log, the bibliography, and the generated result files in `paper/quantum/results/` and `paper/analysis/`. No repository file was edited.

## Summary (3 sentences)

The paper asks how accurate singlet–triplet gaps must be for a golden-rule/MLJ ISC screening decision on BODIPY dyes to be right (about 0.08 eV for a factor-2 rate, 0.016 eV for a near-resonant upper triplet). It then checks a B3LYP/ΔSCF pipeline and (12e,11o) multireference models against three measured rates, implements the golden-rule rate as a QSVT spectral filter, and costs a per-dye screening decision on a fault-tolerant computer (10^11–7×10^15 Toffolis). The work is careful and unusually candid about its limitations, and most numbers check out against the tables and generated outputs. However, the SI decision-cost table is stale and contradicts the main text, one headline number in the Conclusions has no support, a `\todo` still shows in the PDF, and the central message is weakened by two problems the paper itself documents: the spin–orbit uncertainty alone uses up the factor-2 budget, and the costs depend on a window-dependent spin–orbit norm.

## Recommendation

**Return to the author for revision before review.** The blocking items (findings 1, 2, 5) are text or regeneration fixes that take hours, not new science. They should be fixed before referees spend time on the paper, because a referee who compares Table 4 with Table S9 will lose confidence in every number. Findings 3 and 4 concern how the result is framed. Fixing them mostly means rewriting text, plus one inexpensive calculation the existing code can already do. Desk rejection is not warranted: the paper fits JCTC (an accuracy-requirement analysis plus resource estimation for a chemical observable), and the work is substantive and reproducible.

---

## Findings

Legend: **[T]** = text-only fix; **[R]** = regenerate an existing table or figure with existing scripts (no new science); **[N]** = requires new work.

### Major

**1. SI Table S9 (all-dye decision costs) is stale and contradicts main Table 4 and the main-text ranges.** [R]
- Location: si_tables.tex:70–102 (si.pdf p.10, Table S9) vs main.tex:647–656 (Table 4, main.pdf p.22) and main.tex:598–612.
- Evidence: the main Table 4 rows match the current generator output `paper/quantum/results/decision_table.tex` and `decision_cost.md` (both dated Nov 22 16:57). Table S9 was written by `analysis/make_si_tables.py` (si_tables.tex dated before the last decision_cost run) and holds the older values:
  - I2-BDP: D1 MLAE 1.8e14 (main) vs 1.6e14 (S9); canonical 1.6e15 vs 1.4e15; channel 1.1e13 vs 1.2e13; **D2 8.8e13 vs 1.8e13**.
  - I-BDP: D2 9.5e13 vs 1.6e13. BDP: D2 3.1e11 vs 9.6e10. Br2-BDP: D2 9.5e13 vs 2.8e14. aza-BDP-IPh2: D2 6.4e14 vs 4.7e14 (the abstract's "up to 6×10^14" relies on the main value).
  - aza-BDP-BrPh2: D2 4.0e11 in S9, which contradicts main.tex:604 "1×10^12–1×10^14 for the eight heavy-atom dyes outside the zone" (current output: 1.3e12).
  - Light-atom D2 in S9 is 5.1e10–5.3e11; the main text (l.603) says 1×10^11–2×10^12 (current output: 1.3e11–2.2e12).
  - The channel-route speed-ups from S9 are 9–21×, against "12–28 times" in the abstract (l.58) and l.618.
  - The S9 caption says "D2 … by MLAE", but main Table 4 says "best estimator", and sampling is best for aza-BDP-Ph4 and BDP-0 (decision_cost.md).
  - S9 also lists D1 costs (up to 4×10^19) for dyes with k_ISC < 10^6 s^-1, which the generator and main Table 4 exclude by design.
- Fix: rerun `make_si_tables.py` after `decision_cost.py` and rebuild si.pdf. Use the same D2 rule ("best estimator") in both tables. Add a build check (for example a Makefile dependency) so that Table 4 and Table S9 come from the same JSON.

**2. The Conclusions state a factor that nothing in the paper supports: "the rate changes by a factor of 3–40 between spin–orbit operators".** [T]
- Location: main.tex:708 (main.pdf p.24).
- Evidence: the results give a factor of 2.8 for the dominant I2-BDP channel (l.140, "2.8 times more rate"; Limitations l.696 "a factor of 2.8"), a factor of about 5 for the I2-BDP-0 T2 channel (l.141), a factor of about 3 overall (l.147, "uncertain by a factor of about 3"), and a factor of 3 between the 0.6–9×10^8 and 1.8×10^8–2.6×10^9 ranges (l.309–310). No "40" appears anywhere in main.tex, the SI, or the analysis notes (searched). Table S6 gives channel ratios of up to about 90 for the I2-BDP-0 T5 channel, but that channel is explicitly excluded.
- Fix: replace with the supported number, for example "by a factor of about 3 for the dominant channel (up to 5 in a minor channel)". Do not let the Conclusions introduce numbers the Results do not contain.

**3. The central message ("gaps to about 0.08 eV") is not the binding requirement by the paper's own budget, and the 0.08 eV figure comes from a model the paper says is inadequate.** [T + small N]
- Location: abstract l.43–46; Conclusions l.704; Methods l.186–188; Results l.369–377; si_fcwd.tex.
- Evidence:
  - (a) Methods l.187–188 state that "an SOC model good only to a factor of 3 in k_ISC would leave no room for gap errors at F=2". The paper then finds that its own SOC treatment is uncertain by a factor of about 3 (l.147) and differs from SOMF by 68% in V1 (l.376–377). Read together, these say that for this pipeline the coupling, not the gap, is the limiting input. The abstract and Conclusions present only the gap tolerance.
  - (b) The 0.08 eV value is the half-width of a single resolved vibronic line of the single-mode MLJ model (l.371–374: "its factor-2 tolerance is the half-width of that line, sqrt(4λ_s k_BT ln2)=0.084 eV … the 0.08 eV is model dependent"). The model's reorganization energy (0.24 eV) lies outside both computed estimates (0.12 energy-based, 0.19 harmonic; l.164–166, SI §S12 closing sentence "the two estimates bracket the line-shape uncertainty"). A non-specialist reader of the abstract would take 0.08 eV as a property of BODIPY ISC, when it is a property of the chosen line shape.
- Fix:
  - (i) [T] State in the abstract and Conclusions that coupling and gap errors share the same factor-2 budget, and that with the present SOC uncertainty the coupling is the tighter constraint.
  - (ii) [N, small] The multimode Kubo line shape is already implemented (`fcwd_multimode.py`) for BDP and I2-BDP. Recompute the factor-2 gap tolerance with that density of states and report it next to the MLJ value. If the two agree, the headline is strengthened. If not, the abstract should quote the multimode number or the range.

**4. The decision costs depend on a spin–orbit norm that the paper shows is an artefact of the active-space choice, yet the abstract quotes them as results.** [T; optional N]
- Location: abstract l.52–58; l.574–591; Discussion l.683–688; Table 4.
- Evidence:
  - The cost scales as α_SO/V_eff (l.684). α_SO for I2-BDP is 12 006 cm^-1 (canonical window, B3LYP orbitals, Table 4), 8315 cm^-1 (same window, HF orbitals, Table 2) and 299 cm^-1 (AVAS π space, Table 2), a factor of about 40.
  - The canonical window "holds iodine orbitals by accident" (l.590–591) and is shown in Section 3.3 not to be the chromophore π space.
  - The cost inputs also mix the pipeline rate with the CAS norm. With the CASCI energies of the same Hamiltonian, p changes by a factor of 1200 or 27 (l.579–584).
  - So the dye-by-dye Toffoli counts (10^11–7×10^15) belong to a Hamiltonian that the paper itself argues is not the right one. The robust result is the scaling p = γ²(V_eff/α_SO)² and its decomposition into α_SO/‖χ‖ and ‖χ‖/V_eff.
- Fix: [T] Lead the cost section, abstract and Conclusions with the scaling law and the α_SO/V_eff ratio. Present the dye-by-dye Toffoli counts as illustrative for the canonical window, with an explicit uncertainty of at least an order of magnitude from the active-space choice. [N, optional, would strengthen the paper] Cost one dye (I2-BDP) in a controlled active space that contains the SOC-carrying orbitals (AVAS π plus iodine 5p), which is exactly what Discussion l.678–679 recommends.

**5. A leftover `\todo` renders in red in the compiled PDF; the archival DOI is missing.** [T/admin]
- Location: main.tex:737–738; main.pdf p.25: "[TODO: Archive a tagged release on Zenodo and cite the DOI.]"
- Fix: archive the tagged release, cite the DOI in the Data and Software Availability section, and remove the `\todo`. Without an archived version, the GitHub URL alone does not meet JCTC's data-availability expectations.

**6. The choice of benchmark dyes is not explained (7 of the 10 dyes in ref Ly2021 are used).** [T; possibly N]
- Location: main.tex:261–266; bibliography entry 22 (title: "…in a series of **ten** methylated-meso-phenyl-BODIPY dyes").
- Evidence: the benchmark rests on three measured ISC rates, two of them for one π system (l.690–691), and contribution C4 (l.101–104) calls it "an experimental benchmark of ISC rates and triplet yields". The paper uses BDP, I-BDP, I2-BDP, I3-BDP, BDP-0, I-BDP-0 and I2-BDP-0 from Ly et al. and does not say why the other three were left out. If any of them has a measured k_ISC or Φ_T, leaving it out matters given how thin the validation set is.
- Fix: [T] State the selection criterion. [N] If the omitted dyes have measured rates or yields, add them. This would be the single most effective way to strengthen the classical validation.

### Minor

**7. Wrong cross-reference: "Table S9 for all dyes" for α_SO, but Table S9 has no α_SO column.** [R/T] main.tex:575 (main.pdf p.19). Table S9 columns: k_ISC, p, D1 ×3, D2 calls, D2 Toffolis. Fix: add α_SO and V_eff columns to S9 when it is regenerated (finding 1). They are in decision_cost.md.

**8. Table S7 lists unshifted TDA T2/T3, but its caption does not say so, so gaps read from it conflict with the main text.** [T] si_tables.tex:7–9 (si.pdf p.9). From Table S7, T2 − S1 of aza-BDP-Br2 is 1.96 − 1.80 = 0.16 eV, whereas the main text (l.394) and Tables S8/S10 give +0.11 eV (the shifted pipeline value; 0.16 eV is the "no shift" case in Table S10). The same applies to I2-BDP: 0.34 from S7 vs 0.315 in Table S11. Fix: write "raw TDA T2 and T3 (not shifted by the ΔSCF T1 correction used for the rates)", or list the shifted values.

**9. T5 position and V² ratio are misstated in the SOC benchmark paragraph.** [T] main.tex:138 (main.pdf p.4): "T5, 2 eV above S1, is left out; its V² differs by a factor of 1.7". Table S6: T5 lies at 3.80 eV (I2-BDP) and 2.83 eV (I2-BDP-0), which is 1.5 and 1.1 eV above the ΔSCF S1 (2.26, 1.70) and about 2.1 eV above T1 of I2-BDP. The V² ratio is (1119/848)² = 1.7 for I2-BDP but (29.0/3.0)² ≈ 90 for I2-BDP-0. Fix: give the correct energy reference, and state the ratio for each dye.

**10. Main Table 1 omits two of the nine benchmark dyes (BDP-0, I-BDP-0), and Table 4 omits seven dyes without saying why.** [T] main.tex:333–342 has 7 rows, but the text speaks of nine validation dyes (e.g. l.420, l.690). Table 4 (l.648–655) includes BDP and I-BDP (D2 only) but leaves out the core, aza-BDP-Ph4, BDP-0, I-BDP-0, Cl2-BDP, I1-BDP and aza-BDP-BrPh2, although some of these have k_ISC ≥ 10^6 s^-1 and fall under the caption's own D1 criterion. Meanwhile l.574 speaks of "the heavy-atom dyes of Table 4", and BDP is in that table. Fix: give the selection rule in the caption, or include all rows that meet it.

**11. Undefined or overloaded symbols.** [T]
- λ′ (main.tex:228) is defined only in SI §S4.
- k_th first appears as a column header in Table 4 (l.644). "threshold rate" is used at l.219 without the symbol.
- C_SO (l.226) is never defined in the main text or SI. decision_cost.md shows it is the H_SO basis rotation, 4.5×10^3 Toffolis.
- α in SI §S4 ("divided by α/Δ", si.tex:166) is used before its definition in §S5.
- "ECP" (l.146) is never introduced as an abbreviation.
- λ is used for both the THC one-norm (Ha) and the reorganization energy λ_s, often in the same sentence (l.228–229). Rename the one-norm (e.g. Λ or λ_THC).
- δ is used for the gap error and, in "δ ln V" (l.187), for a relative coupling error.

**12. Figure 2 axis label does not match the definition of p.** [R] figures/decision_cost.pdf: x-axis "success probability p = k_ISC/K′"; caption (l.664) and text (l.576) define p = γ²k_ISC/K′. Fix: relabel the axis.

**13. Conclusions and abstract do not match the Results in three places.** [T]
- (a) The Conclusions (l.706–709) cite only the multimode agreement ("reproduces to within a factor of 2") and leave out the energy-based-reorganization result (4–60× below experiment, l.309–310). The Limitations (l.692–693) say the line shape "alone leaves an order of magnitude in k_ISC".
- (b) Section 3.2 calls window membership "the tightest requirement for red-absorbing dyes" (l.386), but this finding appears in neither the abstract nor the Conclusions.
- (c) The abstract's "up to 270 times cheaper than sampling for heavy-atom dyes" comes from I-BDP (decision_cost.md: 2.6e16/9.5e13), a negative with k_ISC ≈ 3 s^-1 whose iodine sits on the meso-phenyl. For the ISC-active heavy-atom dyes the gain is about 4–140×. Fix: state the range for ISC-active dyes, or say which dye gives the 270.

**14. Limitations say "one conformer" but three were computed for the distyryl dyes.** [T] main.tex:694 vs l.350–351 and Table S13. Fix: "one conformer (three for the distyryl dyes, Table S13)".

**15. The SI list in the main text omits SI Section S12 (multimode line shape, Tables S14–S15)**, although the main text cites S12 twice (l.166, l.313). [T] main.tex:725–730.

**16. The SI claims every table and figure names its script; several do not.** [T] si.tex:41 says "every table and figure lists the script that produced it". Tables S1, S5, S7–S9 and Figures S1, S3, S4 do not, and neither do main Tables 1–3 or Figure 1. Fix: add the script names, or soften the claim to "the repository README maps each table and figure to its script".

**17. Citation problems.** [T]
- ref 9 (Zhou2026PRA) is given as "Phys. Rev. A 2026, arXiv:2512.15889" with no volume or article number. Give the full citation if published, otherwise cite it as a preprint and rename the key.
- ref 17 (Brassard 2002) lacks the chapter title ("Quantum amplitude amplification and estimation").
- ref 49 (Chibani 2013) lacks the page range.
- ref 29 (PySCF) uses "et al.".
- l.74: Gómez-Bombarelli 2016 screened TADF OLED emitters, not "such dyes" (BODIPY photosensitizers). Reword to "high-throughput TDDFT screening of organic emitters is established".
- l.680: "downfolded and transcorrelated Hamiltonians" carries no citation.

**18. Typesetting.** [T/R]
- Table S14 is 57 pt too wide (si.log:616, si.pdf p.14).
- Tables S14–S15 (si_fcwd.tex) have no `\label`, use `\hline` instead of booktabs, and print "-0.718" with a hyphen instead of a minus sign.
- SI §S3 is a heading with only a figure; its text is in §S4. Merge the two sections.

**19. No competing-interest statement or funding statement.** [T] ACS requires a conflict-of-interest disclosure. Add "The author declares no competing financial interest" and a funding line, if applicable.

**20. Abstract length and readability.** [T] The abstract is about 297 words (main.tex:40–58), which is long for JCTC (check the current author guidelines). It packs in about 20 numerical claims, and the sentence at l.46–49 ("…within a factor of 2, with an assumed reorganization energy twice … but 1.3 times …, and with an effective spin–orbit operator, and it fails for…") is hard to parse. Contribution C2 (l.98–100) is similarly garbled. Fix: cut the abstract to the four main findings, roughly one sentence each, and move the qualifiers into the body.

### Scope, title, structure (third-referee view, broad theoretical chemist)

**21. Title and framing.** [T] The title promises a general answer ("How accurate must excited states be…"), but the answer depends on the dye family, the line shape and the coupling model, and about half the paper is fault-tolerant resource estimation, which the title does not mention. The running title, "Accuracy budget and quantum cost of photosensitizer screening", describes the content better. Suggestion: "Accuracy Budget and Fault-Tolerant Cost of Screening BODIPY Triplet Photosensitizers". The final sentence of the Conclusions ("a general recipe") is justified as a *recipe*, but the numbers are not general, and the paper should say so explicitly.

**22. Readability for a non-specialist.** [T] The paper is internally consistent: I checked about 60 quoted numbers against Tables 1–4, S1–S15, accuracy_budget.md and decision_cost.md, and apart from findings 1, 2, 8 and 9 they agree. It is nonetheless very dense. Sections 3.5–3.6 use THC, MLAE, DKLR, SPRT, γ², g_S, K′, p_th, indifference zones and Toffoli counts with little explanation. The following would help a non-quantum reader a great deal:
- one schematic figure of the workflow: pipeline → budget → filter → decision cost;
- a sentence on why Toffoli counts are the cost metric;
- a short "key results" box or table listing each headline number with its section and model assumptions.

Section 3.3 (three paragraphs, two tables) could move partly to the SI. Overall length (32 pp. double-spaced, 4 tables, 2 figures, 1 chart) is acceptable for JCTC.

**23. Reproducibility.** [T/admin] Strong overall. Every script named in the text exists in the repository (checked: tests/test_tdwin.py, paper/references/run_*.py, paper/quantum/*.py, paper/analysis/fcwd_multimode.py, validation/run_xc_scan.py, validation/conformers.py, environment.yml), and software versions are given. Remaining gaps: the archive DOI (finding 5); the regeneration chain that produced the stale Table S9 (finding 1); and the empirical MLAE constant (45) together with the assumed γ² = 0.5 and 10^4-CSF trial state. All are disclosed, but they should also appear in the Limitations paragraph, which currently mentions only "an assumed trial-state overlap".

---

## Checks that passed (no action)
- All main-text SI references resolve to the right content in si.pdf: Table S1/Fig S1 (absorption), S2–S4 (canonical window), Fig S3 (amplitude estimation), Fig S4 (MLJ gap), Table S5 (QSVT), Table S6 (SOC), Tables S10–S13, Section S4 (decision model), S5 (QSVT details), S12 (multimode). The exceptions are "Table S9 … α_SO" (finding 7) and the stale S9 values (finding 1).
- No undefined references or citations in main.log or si.log, and no BibTeX warnings. The main-text table, figure and equation numbers that the SI hard-codes (Table 2, eq 1) are correct.
- Spot-checked numbers that agree:
  - abstract π weights 2.2–6.7; NEVPT2 gap overestimate 0.33–0.39 eV;
  - tolerances 0.079–0.083 / 0.062–0.083 / 0.18–0.38 eV (Table S8); k_BT ln2 = 0.018 eV; 0.084 eV line half-width;
  - quadrature 0.032 eV; p ratios 1200 and 27; α_SO/‖χ‖ = 11 and ‖χ‖/V_eff = 270; α_SO/V_eff = 150–1.3×10^4;
  - D1/D2 ranges in main l.598–614 against decision_cost.md; channel route 12–28×; 6 years at 10^6 Toffoli/s; 3–60× the proxy cost;
  - SOMF ratios 2.8× and 5×; class margins (7 of 9 > 0.5 eV, two ≈ 0.4 eV); window margins.
