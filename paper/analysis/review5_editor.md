# Round 5: pre-review read (associate editor) and third report (broad theoretical chemist)

Manuscript: "How Accurate Must Excited States Be to Screen Triplet Photosensitizers?" I read main.tex (762 lines) and main.pdf (33 pp.), plus si.tex with si_tables.tex and si_*.tex and si.pdf (23 pp.). I checked them against `pdftotext` output, main.log and si.log (both built 2026-11-25 19:23), main.bbl, `quantum/results/decision_cost.md`, `analysis/tolerance_multimode.md`, `analysis/fcwd_multimode.{json,md}` and `analysis/fcwd_table.py`. No repository file was changed apart from this report.

## Summary (3 sentences)

The revision fixes the two blocking numerical problems of round 4. Table S9 has been regenerated and now matches Table 4 and every cost range in the text. The unsupported "factor 3–40" is gone. The paper also adds a multimode line-shape tolerance (0.05–0.10 eV) and a vertical-gradient treatment of the aza-BDP-Br2 T2 channel. The new text mostly agrees across abstract, Results and SI, but it brings four new problems. The adiabatic T2–S1 shift of the top candidate (−0.04 eV) is smaller than the model's own error on S1 relaxation (0.07 eV). The Methods still justify vertical gaps for all dyes, which the new aza-BDP-Br2 result contradicts. The "k_BT ln2 whatever the line shape" argument is contradicted by the paper's own Table S16. The Conclusions still describe the line shape more favourably than the Results do. The red `\todo`/missing DOI is still in the PDF, Table S14 is now cut off at the page edge (140 pt overfull), and about a third of the minor round-4 items (symbols, Table 1/4 row selection, Fig. 2 axis, COI, citations, dye selection) were not addressed.

## Recommendation

**Return for one short, text-only revision before review (about a day of work), then send it out without another editorial round.** The science is ready for referees. Only the red TODO in the PDF (new finding 1) and the truncated SI table (finding 2) keep it from going out as is: a referee who sees either will read the rest with less trust. Findings 3–6 are framing fixes that a referee would otherwise raise. None needs new computation, though finding 3 would gain from a small one. Desk rejection is not warranted.

Legend: [T] text only; [R] regenerate with existing scripts; [N] new work.

---

## Status of the round-4 findings

| # | Round-4 finding | Status | Evidence / what remains |
|---|---|---|---|
| 1 | SI Table S9 stale, contradicts Table 4 | **Resolved** (residual caption issue) | All S9 D1/D2/channel values now equal Table 4 and decision_cost.md (e.g. I2-BDP D2 8.8e13, aza-BDP-BrPh2 1.3e12, aza-BDP-IPh2 6.4e14), and the channel speed-ups are 12–28×. Remaining: the S9 caption (si_tables.tex:74) still says "D2 … by MLAE", while Table 4 (main.tex:656) says "best estimator". For aza-BDP-Ph4 and BDP-0 the two differ (MLAE 1.6e11/1.8e11 vs sampling 1.3e11/1.6e11, decision_cost.md D2 table). S9 also still lists D1 costs up to 4.7e19 for dyes with k_ISC < 10^6 s^-1, which decision_cost.md and Table 4 exclude by design. Fix: one caption sentence. |
| 2 | "factor of 3–40 between SOC operators" unsupported | **Resolved** | main.tex:728 now says "a factor of about 3 between spin–orbit treatments", which matches l.143 and l.152. |
| 3 | Gap tolerance not binding; 0.08 eV is model-dependent | **Mostly resolved** | (ii) Multimode tolerance computed and reported (l.384–388, abstract l.44). (i) The Results now say that SOC and gaps share the factor 2 (l.389–391). The abstract does not say this. The Conclusions (l.728–729) say only that the agreement "does not yet validate the coupling" and do not say that, at the present SOC uncertainty, the coupling uses up the budget. See also new finding 5 (the 0.05–0.10 range and its λ_out dependence). |
| 4 | Costs depend on an active-space-dependent α_SO | **Partially resolved** | The abstract now says "a spin–orbit norm that varies 40-fold … so they measure how the decision cost scales rather than a quantum advantage" (l.57–59), and l.606–607 carries the caveat. The Conclusions (l.730–732) still quote "1×10^11–7×10^15 Toffolis" without the 40-fold caveat. The optional [N] costing in an AVAS+I 5p space was not done. |
| 5 | Red `\todo`; no archival DOI | **Not addressed** | main.tex:757–758; main.pdf p.26: "[TODO: Archive a tagged release on Zenodo and cite the DOI.]" |
| 6 | Why 7 of the 10 Ly et al. dyes | **Not addressed** | main.tex:268–271 gives no selection criterion. The bibliography title still reads "…a series of ten…" (refs.bib:16–17). |
| 7 | "Table S9 for all dyes" for α_SO, but S9 has no α_SO | **Resolved, by another route** | l.591 now points to `paper/quantum/results/decision_cost.json`. That works for reproducibility, but a reader of the SI still cannot see α_SO and V_eff for the seven dyes left out of Table 4. Adding two columns to S9 would be better. |
| 8 | Table S7 caption did not say T2/T3 are unshifted | **Resolved** | si_tables.tex:7 "raw TDA, before the shift by the ΔSCF correction of T1 that the rates use". |
| 9 | T5 position and V² ratio misstated | **Resolved** | l.141: "1.1–1.5 eV above S1 … factor of 1.7 for I2-BDP and 90 for I2-BDP-0", which matches Table S6. |
| 10 | Table 1 omits BDP-0/I-BDP-0; Table 4 row selection unexplained | **Not addressed** | Table 1 still has 7 rows (l.345–352), while the text refers to "nine" validation dyes (l.394, l.708). Table 4 still includes BDP and I-BDP but omits Cl2-BDP, I1-BDP and aza-BDP-BrPh2, all of which have k_ISC ≥ 10^6 s^-1 and therefore meet the caption's D1 criterion (S9 lists their D1 costs). It also omits the core, aza-BDP-Ph4, BDP-0 and I-BDP-0. |
| 11 | Undefined/overloaded symbols | **Not addressed** | λ′ (l.233) is defined only in SI l.163. k_th (Table 4 header, l.662) is never defined in the main text. C_SO (l.231) is never defined (decision_cost.md: "H_SO basis rotation 4.5e+03"). "ECP" (l.150) is never introduced; l.117 and l.136 spell out "effective core potential" without the abbreviation. SI l.166 uses α before §S5 defines it. "δ ln V" (l.192) reuses δ. λ now means the THC one-norm (l.262), λ_s, λ′, λ_out, λ(X) and λ_S1T1 (Table S14), and wavelength (Table S1), so the overload is worse than in round 4. |
| 12 | Fig. 2 x-axis "p = k_ISC/K′" vs p = γ²k/K′ | **Not addressed** | `pdftotext figures/decision_cost.pdf` still reads "success probability p = kISC/K′"; caption l.682 says p = γ²k_ISC/K′. |
| 13 | Abstract/Conclusions vs Results | (a) **Partially**; (b) **Not addressed**; (c) **Resolved** | (a) The Conclusions now mention the 1.5× uncertainty in λ but still say the multimode model "reproduces [the line shape] to within a factor of 2", whereas Results l.328–330 say "by up to an order of magnitude in k_ISC" and l.320 says "4–60 times below experiment" (see new finding 4). (b) "Window membership is the tightest requirement" (l.400) still appears in neither the abstract nor the Conclusions. (c) The abstract now says "4–140 … for the ISC-active heavy-atom dyes", and the body gives 4–270 for all heavy-atom dyes. Both check against decision_cost.md (I-BDP 274×; aza-BDP-Br2 3.6× to aza-BDP-IPh2 142×). |
| 14 | "one conformer" | **Resolved** | l.712 "one conformer (three for the distyryl dyes)". |
| 15 | SI list omits S12 | **Resolved, new gap** | l.749 now lists "Section S12, Tables S14 and S15". The new Table S16 (aza-BDP-Br2 S1→T2) is neither listed nor cited anywhere in main.tex. |
| 16 | "every table and figure lists the script" | **Not addressed** | si.tex:41 is unchanged. Tables S1, S5, S7, S8, S9 and Figures S1, S3, S4 still name no script, and neither do main Tables 1–3 and Figure 1. |
| 17 | Citation problems | **Not addressed** | Zhou2026PRA still prints "Phys. Rev. A 2026, arXiv:2512.15889" with no volume or article number (main.bbl:111–112), although refs carry a DOI (related_work.bib:18). Brassard2002 still has no chapter title (bbl:180–181). Chibani2013 has no page range (bbl:459). PySCF still uses "et al." (bbl:289). The wording at l.76 ("TDDFT of such dyes … GomezBombarelli2016", a TADF-emitter paper) is unchanged. l.698 "downfolded and transcorrelated Hamiltonians" still has no citation. |
| 18 | Typesetting | **Not addressed; worse** | Table S14 is now **140.7 pt** too wide (si.log:616; si.pdf p.14: the "leading modes" column is cut off at "cm" and loses its third mode). Tables S14–S16 still use `\hline`, have no `\label`, and print "-0.718" with a hyphen. §S3 is still a heading with only a figure. |
| 19 | No COI/funding statement | **Not addressed** | No "competing" or "funding" in main.tex. |
| 20 | Abstract length/readability | **Worse** | The abstract is now about 340 words (was about 297), with about 25 numbers. The sentence at l.46–49 ("… within a factor of 2, with an assumed reorganization energy twice … but 1.3 times …, and with an effective spin–orbit operator, and it fails …") is unchanged, and C2 (l.100–102) is still garbled. |
| 21 | Title/framing | **Not addressed** | Same title. |
| 22 | Non-specialist readability (workflow figure, key-results box) | **Not addressed** | Still 2 figures + 1 chart; main.pdf grew from 32 to 33 pp. |
| 23 | Reproducibility: MLAE constant, γ², 10^4-CSF trial state in Limitations | **Not addressed** | l.718 still says only "an assumed trial-state overlap". |

---

## New findings

### Blocking for send-out (presentation)

**1. The red `\todo` and the missing archive DOI are still in the PDF.** [T/admin] — round-4 #5 carried over.
- Location: main.tex:757–758; main.pdf p.26.
- Fix: archive a tagged release, put its DOI in the Data and Software Availability section, and delete the `\todo`. This is the one item that keeps the paper from being sent out as is.

**2. Table S14 overflows the page by 141 pt and loses data.** [R]
- Location: si_fcwd.tex:6–14 (generated by analysis/fcwd_table.py); si.log:616; si.pdf p.14.
- Evidence: in the rendered table the last column header is cut to "leading modes, cm", and the third leading mode of each row (e.g. "1287 (0.09)") is missing.
- Fix: in fcwd_table.py, use `\small` with `\setlength{\tabcolsep}{3pt}` and a `p{}` column, or print only two leading modes. Switch to booktabs, add `\label`s, and write `$-0.718$`. Also rename §S12, which still says "Multimode line shape for BDP and I2-BDP" although it now covers aza-BDP-Br2 too.

### Major (scientific framing; text-only unless stated)

**3. The adiabatic T2–S1 shift of the top candidate (−0.04 eV) is smaller than the vertical-gradient model's own error on S1, yet the Results state it as a fact.** [T; optional N]
- Location:
  - main.tex:302–304: "the vertical-gradient model … relaxes it 0.04 eV more than S1, so T2−S1 drops from 0.11 to 0.07 eV between the minima";
  - main.tex:416–418;
  - si_fcwd.tex:36.
- Evidence:
  - For aza-BDP-Br2 the vertical-gradient (TDA) S1 relaxation energy from the S0 minimum is 0.153 eV (si_fcwd.tex:36; fcwd_multimode.json `grad:S1` 0.1531). The ΔSCF optimization of the same dye gives an S1 relaxation of 0.08 eV (main.tex:301). The model therefore overestimates S1 relaxation by 0.07 eV, about twice the 0.038 eV difference it reports and four times the 0.016 eV tolerance of that channel.
  - T1 agrees (0.086 vs 0.09 eV), so the discrepancy is specific to S1, the state the gap is measured from.
  - The same pattern holds for BDP and I2-BDP in Table S14: TDA-gradient λ(S1) is 0.105/0.163 against 0.051/0.108 from the ΔSCF minima, and TDA-gradient λ(T1) is 0.074/0.076 against 0.123/0.127.
  - If one combined the ΔSCF S1 relaxation (0.08) with the gradient T2 relaxation (0.19), T2–S1 would be about 0.00 eV, not 0.07 eV.
  - The SI itself says a vertical-gradient line shape "is uncertain by about a factor of 3 in either direction" (si_fcwd.tex:34), but the main Results do not carry that caveat. Only the Limitations do (l.713).
  - Separately, the vertical gap is taken at the xTB geometry and the relaxation from the B3LYP S0 minimum.
- Consequence: the conclusion that the class holds is probably robust, because a T2 closer to resonance only raises the rate. The stated adiabatic gap and the rate range 3×10^9–1.5×10^11 s^-1, however, are not determined to better than about ±0.07 eV in the gap.
- Fix: [T] State that the adiabatic T2–S1 is 0.0–0.11 eV, with the model uncertainty attached, give the S1 cross-check (0.15 vs 0.08 eV) explicitly in the Results and in SI §S12, and keep the class statement. [N, small] A TDA (or TDDFT) T2 geometry optimization with root following, or at least a TDA S1 optimization to put the S1 error on the same footing, would settle it.

**4. The Methods justify vertical gaps for all dyes on grounds the new result contradicts.** [T]
- Location: main.tex:164–167, "differ by +0.03, −0.01 and +0.02 eV … well inside the tolerance of Section 3.2, so we use vertical gaps for all dyes".
- Evidence: these are S1–T1 shifts. For the one dye where the rate runs through a near-resonant upper triplet, the paper's own estimate of the vertical-to-adiabatic change of the gap that matters is −0.04 eV (l.303–304). That is 2.5 times the 0.016 eV tolerance of that channel (l.410). The BDP S1–T1 shift (+0.03) also exceeds 0.016, and l.299–300 acknowledges this.
- Fix: rewrite as "… within the 0.08 eV tolerance of the T1 channels; for the near-resonant T2 channel of aza-BDP-Br2 the vertical-gradient estimate is −0.04 eV, beyond the 0.016 eV tolerance (Section 3.2)."

**5. "For an uphill channel the rate falls as e^{−ΔE/k_BT}, so this tolerance approaches k_BT ln2 = 0.018 eV whatever the line shape" is not true in general, and the paper's own Table S16 contradicts it.** [T]
- Location: main.tex:410–412; the abstract (l.45, "close to k_BT ln2") depends on the same argument.
- Evidence:
  - In a Gaussian (Marcus) line, d ln k/dΔE = −(ΔE+λ)/(2λk_BT). This is steeper than 1/k_BT whenever ΔE > λ.
  - The single-mode model happens to have ΔE ≈ λ_s (0.11 vs 0.10 eV), which is why it gives about 0.016 eV.
  - With aza-BDP-Br2's own S1/T2 modes the mutual reorganization energy is 0.014 eV (si_fcwd.tex:36). Table S16 gives rate ratios between ΔE = +0.069 and +0.107 eV of 7.0, 6.8 and 3.8 for λ_out = 0, 0.05 and 0.10 eV. These correspond to factor-2 tolerances of 0.038·ln2/ln(ratio) = 0.014, 0.014 and 0.020 eV.
  - The headline number (about 0.016 eV) therefore survives, but for a different reason, and it is not a line-shape-independent bound.
  - Also, tolerance_multimode.md computes the aza-BDP-Br2 tolerance with I2-BDP modes ("aza-BDP-Br2 (I2-BDP modes)", main l.387), although the aza-BDP-Br2 Hessian and vertical gradients now exist.
- Fix: [T] Replace the sentence with: "for an uphill channel the rate falls at least as fast as e^{−ΔE/k_BT} once ΔE exceeds the reorganization energy, so the tolerance is at most of order k_BT ln2 (0.014–0.020 eV with the aza-BDP-Br2 modes, Table S16)". [R] Rerun tolerance_multimode.py with the aza-BDP-Br2 modes.

**6. The Conclusions present the line shape more favourably than the Results do, and the abstract's multimode range covers I2-BDP only.** [T]
- Location: Conclusions l.726–728; abstract l.44; Results l.318–330 and l.384–386.
- Evidence:
  - (a) Conclusions: "a line shape that a parameter-free multimode model reproduces to within a factor of 2 in the rate". This is true only for the ΔSCF-minima variant (ρ_multi/ρ_MLJ = 0.60–1.65, Table S15). The TDA-gradient variant gives 0.21, a factor of 4.8. The model is not parameter-free, since λ_out is varied from 0 to 0.10 eV. The Results say the line shape leaves "up to an order of magnitude in k_ISC" (l.329), and with the energy-based reorganization the rate is "4–60 times below experiment" (l.320).
  - (b) The abstract's "0.05–0.10 eV with a multimode line shape" is the I2-BDP range (0.046/0.079/0.098). BDP gives 0.036–0.069 (l.385–386; tolerance_multimode.md).
  - (c) l.386 says "the 0.08 eV holds in solution and halves in the gas phase". The measurements are in toluene, a non-polar solvent for which λ_out of 0.05–0.10 eV for a ππ* state needs justification. The best multimode agreement with the measured I2-BDP rate is at λ_out = 0 (4.1×10^9 vs 3.5×10^9 s^-1; 8.0×10^9 at 0.05). If the multimode model is taken as the validated line shape, the self-consistent tolerance is about 0.05 eV, not 0.08 eV.
  - (d) The Conclusions give neither the multimode range nor the window-membership result (round-4 #13b).
- Fix:
  - Conclusions: "…with a single-mode line shape whose rate agrees with a multimode model built from computed modes within a factor of 2 (5 with vertical gradients); the line shape alone leaves an order of magnitude in k_ISC";
  - abstract: "0.04–0.10 eV for BDP and I2-BDP, depending on the outer-sphere reorganization";
  - Results: add one sentence justifying the λ_out range for toluene, or state that at λ_out = 0, which best matches the measured rate, the tolerance is 0.05 eV.

### Minor

**7. Table S16 is never cited in the main text and is missing from the SI list.** [T] main.tex:418 cites only "SI Section S12", and l.749 lists "Tables S14 and S15". Fix: cite "Table S16" at l.418 and add it to l.749.

**8. The aza-BDP-Br2 summary numbers are inconsistent.** [T] main.tex:419–420 says "Across these treatments the class … holds (Φ_T = 0.80–1.00) and only its rate is fragile, by a factor of about 50". The treatments in that paragraph include the unshifted T2 (1.4×10^9 s^-1, Φ_T = 0.80, Table S10), so the rate spread is 1.4×10^9–1.5×10^11, a factor of about 100. The 50 refers only to the multimode adiabatic range (si_fcwd.tex:47, where Φ_T = 0.90–1.00). Fix: "by a factor of about 100 (50 within the multimode model)".

**9. The SI multimode text uses "from the minima" for two different quantities.** [T] si_fcwd.tex:34 first quotes the energy-based reorganization at the optimized minima ("0.09 and 0.12 eV") and then, in the same paragraph, "(0.13 eV from the minima)" and "(0.19 eV)" for the harmonic projection of the same minima. Fix: write "harmonic projection of the ΔSCF minima" in the second case.

**10. Stored outputs no longer contain the numbers the paper cites.** [R] `analysis/fcwd_multimode.json` and `.md`, which main l.324 cites through `fcwd_multimode.py`, now hold only aza-BDP-Br2. The BDP and I2-BDP values behind l.170, l.324–327 and Tables S14–S15 are produced only by re-running `fcwd_table.py`. `fcwd_table.py:35,38` hard-codes ΔE = 0.107 eV and V2 = 46.77 cm^-1 instead of reading them from the pipeline run. Fix: rerun `fcwd_multimode.py BDP I2-BDP azaBDP_Br2` so the committed outputs cover every number cited, and read ΔE and V from `screening/runs/...`. Also, `fcwd_multimode.json` reports harmonic relaxation energies of 0.40 eV (S1) and 0.29 eV (T1) for aza-BDP-Br2's ΔSCF minima (`min:S1`, `min:T1`), 4–5 times the energy-based 0.08 and 0.09 eV. The SI should say that the minima-projection model is unreliable for this dye (likely phenyl torsions), which is why only the gradient model is used for it.

**11. Small numerical slip.** [T] l.147: "every I2-BDP-0 channel 7–8% below SOMF". Table S6 gives V2 = 3.2 against 3.5 (8.6%) and V1 = 0.4 against 0.5 (20%, at rounding precision). Fix: say "7–9% for the channels above 1 cm^-1".

**12. Abstract and Conclusions drift.** [T] The abstract (l.57–59) now carries the 40-fold α_SO caveat, but the Conclusions (l.730–732) do not (round-4 #4). Add half a sentence: "…for the canonical window, whose spin–orbit norm, and with it the cost, varies 40-fold with the choice of active orbitals".

### Third-referee view (broad theoretical chemist)

**13. The central message is clear in the body but buried in the abstract.** [T] The Results now support a coherent three-part message:
- (i) a factor-2 rate needs gaps to about 0.05–0.10 eV, and to about 0.016 eV through a near-resonant upper triplet, with coupling and line shape sharing the same budget;
- (ii) no method tested here reaches that accuracy across the set (B3LYP/ΔSCF agrees for one π system, perhaps by compensation; NEVPT2 misses by about 0.35 eV; distyryl ISC is unexplained);
- (iii) the decision cost scales as α_SO/V_eff, so the spin–orbit encoding and channel resolution matter more than walk cost.

Each revision round has added qualifiers. The abstract has grown to about 340 words with about 25 numbers, and the paper to 33 + 23 pp. A non-specialist reading the abstract cannot tell which numbers are results and which are caveats.

Recommendations:
- (a) Cut the abstract to about 200–250 words, one sentence per message above (check the current JCTC limit).
- (b) Add the "key results" table proposed in round 4 (number, section, model assumption). This would let the caveats leave the abstract without being lost.
- (c) Move §3.3's canonical-window paragraph (l.526–535) and the qubit-loading discussion (l.566–574) to the SI.

**14. The physics of the near-resonance argument deserves one sentence of context.** [T] A reader from photophysics will ask why spin–vibronic (Herzberg–Teller and nonadiabatic) coupling between T2 and T1 is not considered for aza-BDP-Br2. With T2 only about 0.0–0.1 eV above S1 and V2 = 47 cm^-1, it is the standard mechanism for such systems. The paper mentions it only for I2-BDP-0 (l.367) and in a Table 1 footnote. Fix: one sentence in §3.2 or the Limitations saying that the aza-BDP-Br2 rate (and the 0.016 eV tolerance) assumes direct SOC with no T2–T1 vibronic mixing.

**15. Reproducibility.** [T/admin] Apart from findings 1 and 10 and round-4 #16 and #23, reproducibility is good. Every script named in the new text exists (fcwd_multimode.py, fcwd_table.py, tolerance_multimode.py, run_adiabatic.py), the SI tables are generated, and the S9/Table 4 chain is now consistent.

---

## Checks that passed (no action)

- **Cross-references** in main.pdf all resolve to the right SI content: Table S1/Fig. S1 (absorption); Tables S2–S4/Fig. S2 (canonical window); Fig. S3 and §S4 (amplitude estimation, decision model); Fig. S4 (MLJ); §S5/Table S5 (QSVT); Table S6 (SOC); S7–S9 (all dyes); S10 (upper triplets); S11 (TDA window); S12–S13 (functional, conformers); §S12 (multimode); "Table 2 of the main text" and "eq 1" in the SI. There are no undefined references or citations in either log, and main.log has no overfull boxes.
- **Numbers checked against tables and generated outputs, all consistent:**
  - Table 4 = S9 = decision_cost.md;
  - D1 7×10^12–9×10^14 (MLAE) and 4×10^13–7×10^15 (canonical);
  - calls 5×10^4–4×10^6; D2 heavy-atom 1×10^12–1×10^14 (eight dyes);
  - sampling/MLAE 4–270 (all heavy-atom) and 4–140 (ISC-active); canonical/sampling 1.2–10;
  - unknown-margin penalty 2–5×; channel route 12–28×;
  - about 6 years and 0.4–50 days at 10^6 Toffoli/s; proxy ratios 10^2–10^5;
  - p ratios 1200/27/800; α_SO/V_eff 150–1.3×10^4; α_SO 40-fold (12006 vs 299);
  - γ² measurement at most 0.4% of D1; qubits 200/220;
  - SOMF 2.8× and 5×, 5.6× measured; V² ratios 1.7/90; T5 1.1–1.5 eV;
  - tolerances 0.079–0.083/0.062–0.083/0.18–0.38 eV; factors 21/77;
  - multimode tolerances 0.046/0.079/0.098 and 0.036/0.069/0.068 and 0.017–0.020;
  - multimode rates 4.1×10^9 and 0.8–1.1×10^10; reorganization energies 0.13/0.19 and Σ S_j 1.15;
  - aza-BDP-Br2 relaxations 0.153/0.086/0.191 and Table S16;
  - window margins 0.00–0.07/0.08/0.15–0.21 eV;
  - AVAS and canonical NEVPT2 error ranges (Tables 3, S2–S4); distyryl T2 positions; conformer factors 47/3; TDA-window shifts −0.010 to +0.004 eV.
