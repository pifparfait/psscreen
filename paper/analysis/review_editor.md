# Handling-editor and generalist referee report

Manuscript: "How Accurate Must Excited States Be to Screen Triplet Photosensitizers? An Accuracy Budget, a
Gap-Resolved Intersystem-Crossing Observable, and Decision-Level Quantum Resource Estimates"
(`paper/main.tex`, SI `paper/si.tex`), target J. Chem. Theory Comput.
Reviewed 2026-11-18 from a fresh build (`/tmp/rev-editor/main.pdf`, 33 pp.; `si.pdf`). Both build with no
undefined references or citations. Line numbers refer to `main.tex` unless marked `si.tex`.

---

## 1. Recommendation

**Major revision before submission (do not submit in the present state).** The paper fits JCTC, and once
the items in Section 2 are done I would send it out for review without hesitation. It is not ready now
for three reasons:

1. There are 10 open `\todo` items. Several of them hold up headline claims: the SOC model is not yet
   benchmarked, NEVPT2 is still running for the π-space dyes, the experimental values are not yet
   re-checked, and the THC error is missing from the budget.
2. The π-space SC-NEVPT2 result (a 0.54–0.59 eV S1 error) is the main support for "dynamic correlation at
   second order is not enough". It contradicts the CASPT2/CAS(12,11) literature that the paper and
   Zhou et al. cite (Momeni & Brown 2015), and the paper does not discuss that conflict.
3. The take-home message is restated four times (abstract, Section 3.3, Discussion, Conclusions). In two
   places it compares an absolute S1 error with a budget that is defined for S–T gaps.

### Fit and novelty

- **Fit.** The scope is right for JCTC: method assessment, an observable design, and resource
  estimation, all anchored to experiment. A JCTC reader will find the accuracy budget (C1) the most
  useful and most novel part. The "decision-linked tolerance" framing is new. Classical sensitivity
  analyses of ISC rates exist (Shi 2022, Marian), and the introduction should name the closest one in one
  clause so a referee does not raise it.
- **Is there a genuine quantum-computing contribution?** Yes, but it is modest, and the paper says so
  honestly (l. 122–123: "We do not propose a new quantum algorithm"). The contribution has three parts:
  - The golden-rule ISC rate as a QSVT spectral filter over all triplets and spin sectors is implemented
    and verified at circuit level (Table 5).
  - The filter is costed like for like against the gap-free proxy with the same THC estimator and
    precision rules. The useful observation is that qubit counts depend on how the rotation angles are
    loaded (l. 541–546).
  - Decision-level costing uses an empirically calibrated MLAE constant, and it shows a real
    algorithmic consequence: p ∝ (‖χ‖/α_SO)², so amplitude estimation is mandatory even for clear
    negatives.

  This is enough for JCTC. It would not be enough for PRX Quantum or Quantum.
- **Positioning against the four key works.** Mostly fair. Fix the following:
  - *Zhou et al. (PRA 2026).* l. 80–85 describe the proxy correctly, but they omit that Zhou et al.
    themselves state the proxy "is not an absolute ISC rate", invoke the golden rule, and intend
    calibration against a reference molecule (`analysis/xanadu_comparison.md`). Add half a sentence, for
    example after "contains the SOC but no energy gap": "; the authors note that it is not an absolute
    rate and is meant to be calibrated against a reference dye". Without it, l. 596–597 ("cannot take
    decision D2 without calibration") reads as a criticism they had already anticipated.
  - *Which (12e,11o) did Zhou et al. use?* The abstract (l. 48–50) and Conclusions (l. 687) say that
    "the (12e,11o) window of canonical orbitals used in earlier resource estimates" is not the π space.
    Zhou et al. took the CAS(12,11) of Momeni & Brown, which is a π active space by construction. Unless
    you have checked that their integrals came from a canonical window, the sentence attributes our
    construction to them. Suggested abstract wording: "A (12e,11o) window of canonical orbitals, the
    simplest way to build an active space of the size used in earlier resource estimates, holds only
    2.5–6.7 of the 11 dipyrrin π orbitals".
  - *Eklund et al. 2026.* The handling is good: l. 94, l. 126 ("as Eklund et al. made for photochemical
    yields") and l. 598–600 ("overlap the lower end of the range"). Their numbers are 10¹³–10¹⁶ Toffolis
    on 3000–8000 qubits. Our minimum is 5×10¹⁰, so "overlap the lower end and extend below it" is more
    accurate.
  - *Zak 2026.* l. 85–86 is accurate and sufficient. Optionally add "without implementation or numerical
    results" so the C2 claim of "implemented, verified and costed" has a contrast.
  - *Nykänen 2024.* l. 87 and l. 128 are accurate. Optionally add one clause: their ΔUCCSD(8,8) MAE of
    0.13 eV for S1, against our 0.055 eV for the classical pipeline. That supports the "classical baseline
    first" argument.

---

## 2. Required changes, by priority

### P1. Must be done before submission

**R1. Reconcile the π-space NEVPT2 result with Momeni & Brown CASPT2 (l. 454–465, Table 4).**
Zhou et al. say that CASPT2 with this active space "agrees closely with experiment". We report SC-NEVPT2
in a π (12e,11o) space 0.54–0.59 eV too high. A referee will ask why. The likely reasons are:

- the def2-SVP basis (excitation energies fall with basis size; Momeni & Brown used larger bases);
- SC-NEVPT2 versus (IPEA-shifted) CASPT2;
- state averaging;
- geometry.

Required:

- cite Momeni & Brown's CAS(12,11)-CASPT2 numbers explicitly in this paragraph;
- run at least one basis-set check: BDP SC-NEVPT2 with def2-TZVP, or with cc-pVTZ on the AVAS orbitals;
- state whether the 0.33–0.37 eV S1–T1 gap error survives the larger basis.

The gap error matters more than the S1 error for the budget. Table 4 also shows that CASCI already
overestimates the S1–T1 gap by about the same amount as NEVPT2 (1.00–1.10 eV against 1.05–1.11 eV), so
NEVPT2 does not improve the gap at all. Say this. It strengthens the message, but it also points to a
reference or basis problem rather than a lack of dynamic correlation. If the basis check closes most of
the gap, the claim "second-order dynamic correlation in (12e,11o) does not reach the budget" (l. 462–464,
l. 499, l. 643–644, l. 688–689 and the abstract) must be removed or weakened throughout.

**R2. Do not compare S1 errors with a gap budget (l. 487, l. 642, l. 690).** Examples are "CASCI
overestimates S1 by 1.4–2.0 eV, 20–30 times the budget" and "bare active-space energies, which miss S1 by
20–30 times the budget". The budget (Section 3.2) is a tolerance on S1–Tn gaps, so an absolute S1 error
does not enter it. There is also an inconsistency: in the AVAS π space with SA-CASSCF orbitals (Table 4),
bare CASCI misses S1 by only 0.40–0.56 eV, not 1–2 eV. The 1.4–2.0 eV figure holds only for the
canonical window on RHF orbitals. Suggested Conclusions wording (l. 689–690):

> "A gap-resolved ISC observable, implemented with QSVT and verified against the exact golden rule,
> ranks dyes wrongly on bare small-active-space energies, whose S1–T1 gaps are 1.7–2.2 eV too large."

The 1.7–2.2 eV follows from Table 5: 2.3–2.8 eV against about 0.6 eV. In the Discussion (l. 642), replace
"bare active-space Hamiltonians miss S1 by 1–2 eV" with "bare active-space Hamiltonians overestimate the
S1–T1 gaps by 0.3–0.4 eV even in a controlled (12e,11o) π space, and by far more in small or canonical
windows".

**R3. Close or explicitly scope every TODO.** See Section 4. The three that affect claims are:

- l. 152 (SOC benchmark): every rate rests on an untested SOC model;
- l. 268 (experimental values): the benchmark itself;
- l. 464 (π-space NEVPT2 still running): Table 4 has an empty I3-BDP row.

**R4. Fix the Table 5 caption–content mismatch (l. 521–535).** The caption lists three rates: the QSVT
circuit, "the same polynomial on the exact spectrum", and the exact golden rule. The table has two rate
columns. The text (l. 505–507) claims the circuit matches the polynomial "to all printed digits" and the
exact rate to 0.03 % (2 % for the smallest). The polynomial-on-exact-spectrum column is therefore missing,
and the reader cannot check the first claim. Either add the column or change the caption to "from the QSVT
filter circuit and from the exact golden rule" and give the circuit-versus-polynomial agreement as a
number in the text (for example "to 10⁻⁹ relative").

**R5. Explain the two α_SO values for the same dye.** Table 3 (canonical window, "Hartree–Fock orbitals")
and Table 6 / `quantum/results/alpha_so.json` disagree:

| Dye | Table 3 (cm⁻¹) | Table 6 / alpha_so.json (cm⁻¹) |
|---|---|---|
| BODIPY core | 133 | 199 |
| BDP | 27 | 39 |
| I-BDP | 5579 | 6553 |
| I2-BDP | 8315 | 12006 |
| I3-BDP | 11421 | 12372 |

l. 564 hints at the cause ("with the same Hartree–Fock orbitals, 8300 cm⁻¹"), which suggests Table 6 uses
Kohn–Sham orbitals. Neither the caption of Table 6 nor Section 2.4 says which orbitals were used. Add
"α_SO from B3LYP orbitals" (or whatever is correct) to the Table 6 caption and to l. 240–241, and state
in l. 563–565 that the comparison uses HF orbitals for both spaces.

**R6. Number figures and tables in order of first citation (a JCTC requirement).**

- Figure 2 (`fig:tolcurve`) is first cited at l. 183, before Figure 1 (`fig:mlj`, first cited at l. 371).
- Table 2 (`tab:isc`) is first cited at l. 182, before Table 1 (`tab:abs`, first cited at l. 280).

Either drop the forward references from the Methods (l. 182–183: "...within that factor (Section 3.1),
and report the tolerance as a function of F and λ_s in Section 3.2"), or swap the float order.

**R7. Add a Chart of structures.** The paper names 15 dyes with ad hoc labels:

- underscores from the code (`azaBDP_Br2`, `azaBDP_IPh2`, l. 390, l. 394, l. 591, Table 6);
- confusable pairs (I-BDP = p-iodophenyl and I1-BDP = mono-iodo core; aza-BDP-Ph₄ versus azaBDP\_IPh2).

None of the six screening candidates is defined structurally anywhere. Add Chart 1 with all 15 structures
and one consistent label set (for example aza-BDP-Br₂ and aza-BDP-(IPh)₂, without underscores), and use
those labels in all text, tables and the SI.

**R8. Resolve the budget inconsistency at l. 369–373.** The text says "a 0.1 eV error changes k_ISC by a
factor of 1.7–2.7 for all nine dyes" and also that "keeping k_ISC within a factor of 2 requires errors
below 0.062–0.083 eV". If a dye changes by only 1.7× at 0.1 eV, its factor-2 tolerance must exceed 0.1 eV
unless the other sign is worse. Say "for the more sensitive sign" or give the 1.7–2.7 range for both
signs. The two sentences as written cannot both be literally true.

**R9. Take-home message: state it once.** Abstract l. 57–59, Discussion l. 644–646 and Conclusions
l. 694–697 carry nearly identical "useful quantum target" sentences, and the Conclusions (288 words)
re-list almost every number in the abstract. Proposed fix:

- Keep the full statement in the abstract and in the opening sentence of the Discussion.
- Cut the Conclusions to about 120 words that say what changes for the field. Suggested text:

> "Screening triplet photosensitizers requires singlet–triplet gaps to about 0.07 eV, an accuracy a
> TDDFT/ΔSCF pipeline reaches for compact BODIPYs but neither a bare nor a second-order-corrected
> (12e,11o) active space does. A gap-resolved, golden-rule ISC observable can be built and verified
> with QSVT at a per-call cost comparable to gap-free proxies, but a screening decision costs 10³–10⁶
> calls because the rate is a small fraction of the spin–orbit block-encoding norm; amplitude
> estimation is therefore mandatory, and a tailored spin–orbit encoding is the main cost lever. The
> accuracy budget, perturbing the rate model until the decision changes, is a general recipe for
> setting the target of any rate-based quantum screening."

**R10. Abstract length and density.** The abstract is 263 words, which is likely over the JCTC limit
(ACS guidance is 250 words; please check the current author guidelines). It carries 14 numbers. Cut it
to about 200 words. Suggested cuts:

- the "(12e,11o) window ... 2.5–6.7 of the 11 dipyrrin π orbitals" sentence, which can move to the
  results;
- "with tensor-hypercontraction qubitization";
- "about 2 mHa per state".

Define TDDFT and ΔSCF at first use, or rephrase as "a time-dependent DFT pipeline with ΔSCF-corrected
lowest states".

### P2. Should be done

**R11. Title.** The title is 27 words in three clauses. Suggested options:

- "How Accurate Must Excited States Be to Screen Triplet Photosensitizers? Implications for Quantum
  Algorithms" (15 words);
- "An Accuracy Budget for Screening Triplet Photosensitizers and the Cost of Deciding on a Quantum
  Computer".

The short title (l. 25) is fine.

**R12. The mis-ranking claim needs its scope (abstract l. 52–53; l. 510–518).** "On bare active-space
energies it ranks the dyes in the wrong order" rests on two dyes, the BODIPY core and I2-BDP, at CAS(6,6),
where S1–T1 is 2.3–2.8 eV too large. Two options:

- Scope the claim: "in small active spaces".
- Better, evaluate the exact golden rule classically (no circuit needed) on the AVAS (12e,11o) CASCI and
  NEVPT2 energies already in Table 4 and SI Table S2, for all dyes that have them. This is cheap and
  would turn a two-dye anecdote into a table column.

**R13. Window-membership paragraph (l. 383–387).** "Their computed S1 lies 0.15–0.21 eV inside, so a
sharp window classifies all three wrongly": "three" refers to the distyryl dyes, but the preceding
sentence also mentions aza-BDP-Ph₄. That dye's computed S1 (678 nm, 1.83 eV) is 0.08 eV inside the edge,
not 0.15. Make the sentence explicit, for example: "The computed S1 of the three distyryl dyes lies
0.15–0.21 eV inside the window, so a sharp window classifies all three wrongly; aza-BDP-Ph₄ (0.08 eV
inside) is classified correctly only because it sits on the edge."

**R14. l. 559–561 cites a dye that is not in the table.** "I-BDP and I-BDP-0, which have no ISC, reach
700 (Table 6)": I-BDP-0 is not in Table 6. Its ratio of 2673/4 ≈ 670 is in `decision_cost.md` only.
Either add the row or drop I-BDP-0 from the sentence. The Table 6 caption should also say that it shows
8 of the 15 dyes and that the full table is in the SI. Without the full table, the light-atom D2 range
"5×10¹⁰–5×10¹¹" (l. 590) and the D1 range (l. 583) cannot be checked from the paper: of the light-atom
negatives, only BDP appears.

**R15. Number of calls (abstract l. 55; Conclusions l. 693).** The text gives "10³–10⁶ calls", but the
range in `decision_cost.md` is 4.3×10² (aza-BDP-Ph₄, BDP-0, D2) to 2.9×10⁶ (I2-BDP-0, D1). Use "several
hundred to a few million calls" or "10²·⁶–10⁶·⁵". Also make the per-call cost uniform: the abstract and
l. 692 say "about 2×10⁸", while l. 550 and Table 6 say "1.4–2.1×10⁸". Both are fine, but use
"(1.4–2.1)×10⁸" at first mention.

**R16. Qubit count: one number.** The Conclusions give "131–374 logical qubits", Table 6 and l. 550 give
158, and the abstract gives none. Suggested Conclusions wording: "158 logical qubits (131–374 depending on
how rotation angles are loaded)".

**R17. Forward reference to results "below" (l. 462–463).** "...confirms the conclusion of the
canonical-window results below". Either put the canonical-window paragraph first (it is the motivating
observation, and the AVAS run is the controlled follow-up), or move it to the SI (see C3) and rephrase as
"confirms, with more systematic errors, the canonical-window results in SI Tables S1–S3".

**R18. Undefined acronyms.** Expand at first use in the main text (and in the abstract where the term
appears there):

| Acronym | First use | Note |
|---|---|---|
| TDDFT | l. 47 | |
| ΔSCF | l. 47 | |
| QSVT | l. 52, l. 81 | quantum singular value transformation; never expanded |
| SA-CASSCF | l. 50 | |
| SC-NEVPT2 | l. 50 | expanded only at l. 248 |
| CASCI | l. 247 | |
| AVAS | l. 255 | |
| RHF | l. 254 | |
| MINAO | l. 254 | |
| THC | abstract l. 54 | expanded in the abstract, fine |
| TADF | l. 664 | the phrase appears at l. 70 but the acronym is never introduced |
| QPE | l. 653 | |
| DKLR | none | "Dagum–Karp–Luby–Ross" written out, fine |

Also remove QAE from `\abbreviations` (l. 29): it is never used in the text.

**R19. Data availability.** Release the code and data with a tagged version and a DOI (TODO at l. 721).
The repository URL must be public before submission.

**R20. TOC graphic.** JCTC requires one (TODO at l. 36).

---

## 3. Structure and length: suggested cuts

The main text is about 5,900 words of prose by detex (about 7,100 including captions and tables), with
6 tables and 3 figures. JCTC has no hard limit, but the paper argues four things in turn, and every
section repeats the previous conclusion before adding its own. Target: about 5,000 prose words, 4–5
tables, 3 figures plus a Chart.

| # | Cut / merge / move | Lines | Saves |
|---|---|---|---|
| C1 | Conclusions: shorten to about 120 words (R9). Removes the third repetition of all headline numbers. | 683–698 | ~170 w |
| C2 | Discussion, "What the useful quantum target is": drop the first sentence, which recaps Sections 3.1–3.3 (l. 640–644), and start at "We argue that…". Keep the two qualifications. | 640–646 | ~70 w |
| C3 | "Results on the canonical window": move to the SI next to Tables S1–S3. Keep two sentences in the main text: (i) on RHF orbitals CASCI overestimates S1 by 1.4–2.0 eV and NEVPT2 corrects it to 0.06–0.27 eV, but this close agreement is an artefact of the window's mixed character; (ii) for aza-BDP-Ph₄, NEVPT2 puts T2 0.17 eV above S1 against 0.29 eV in the pipeline. Also resolves R17. | 485–499 | ~170 w |
| C4 | "The distyryl failure": keep the finding, the ΔSCF/T2 mechanism and the "three of four methods agree" sentence. Move the PBE0 numbers, the conformer test and the cis/trans SMILES anecdote to the SI. Collapse the list of untested hypotheses to one clause. | 344–364 | ~150 w |
| C5 | Methods, decision-level cost model: move the composition of C_S1 (10⁴ determinants, eigenstate filter, overlap 0.5), the degree fit and the λ′ definition to an SI section "Cost model". Keep the equation C = C_S1 + C_SO + d·C_walk with the values. | 223–232 | ~110 w |
| C6 | "Amplitude estimation is required": move the details of the canonical-QAE grid, the 100-shots-per-depth saturation and the ensemble procedure to the SI with Figure S3. Keep the law check, the two empirical constants and the 5 % validity limit. | 570–581 | ~80 w |
| C7 | Introduction, paragraph 2 (l. 78–101, about 300 w): the classical-side sentence (l. 94–101) duplicates the citations used later at l. 141–142 and l. 284–285. Trim to one sentence. Also merge "three questions" (l. 103–119) with "contributions C1–C4" (l. 123–129): there are currently three questions and four contributions, and the C4 benchmark answers none of the questions. | 94–129 | ~120 w |
| C8 | Merge Figures 1 and 2 into one figure: (a) MLJ curve with channels, (b) tolerance against factor, (c) tolerance against λ_s. Saves a float and fixes R6 at the same time. | 399–415 | 1 figure |
| C9 | Merge Tables 3 and 4 (window character and π-space energies), or move Table 3's AVAS columns into the Table 4 caption/text. Both describe the same (12e,11o) spaces. | 432–483 | 1 table |
| C10 | Table 1: move to the SI with Figure S1 and keep the summary statistics in the text (−0.055 ± 0.02 eV compact; −0.21 to −0.23 eV distyryl). Alternatively merge Tables 1 and 2 into one classical-baseline table: S1 calc/exp, ΔE_ST, V1, k_ISC, Φ_T. Merging is preferable because it removes a table with no information lost. | 291–342 | 1 table, ~0 w |

Estimated saving: about 870 words of prose, 1 figure and 2 tables. Add 1 Chart.

**What should come back from the SI.** Figure S3 (amplitude-estimation simulation) is the only direct
evidence for the statement "amplitude estimation is required", which is in the abstract and the
conclusions. Make it panel (b) of Figure 3 (Figure 3 then shows cost against p next to the error
against calls on the real circuit). Nothing else needs to come back. Table S2 and Table S3 support
side remarks only.

**What the SI is missing.** It should hold:

- the full Table 6 (15 dyes);
- per-dye accuracy-budget tables (the tolerances behind the ranges at l. 372–375);
- the full Table 2 with BDP-0 and I-BDP-0, which are among the "nine validation dyes" but are absent
  from Table 2;
- geometries (xyz);
- the cost-model section (C5);
- the canonical-window narrative (C3);
- the distyryl details (C4);
- the SOC benchmark (R3).

The `suppinfo` paragraph (l. 706–710) must list exactly what is there.

---

## 4. Open TODOs and pre-submission checklist

| Line | TODO | Effect if left open | Action |
|---|---|---|---|
| 36 | TOC graphic | Submission blocked | Make it (tolerance curve plus decision-cost curve, as planned). |
| 152 | SOMF/X2C SOC benchmark of V1 (I2-BDP), V2 (I2-BDP-0) | All rates, the budget's fixed-SOC assumption and α_SO rest on an untested SOC operator. A referee will ask for this first. | **Required.** At least I2-BDP V1 and I2-BDP-0 V2 with SOMF(1X) in ORCA or X2C in PySCF; report the ratio. |
| 268 | Re-check experimental values against Ly 2021 Table 2; aza-BDP-Ph₄ ±5 nm | Benchmark integrity | **Required.** Quick; fix the text "Table 2 of ref 70", which reads ambiguously next to the paper's own Table 2. |
| 287 | Vertical→band-maximum and vertical→adiabatic corrections | The paper admits that the pipeline's agreement may be error cancellation. A referee may push. | Strongly recommended for BDP and I2-BDP. Otherwise reword as a stated limitation and delete the TODO. |
| 464 | I3-BDP π-space NEVPT2, aza-BDP-Ph₄, distyryl | Table 4 has a blank row. The "five dyes so far" in the Limitations (l. 675) is inconsistent with only four NEVPT2 rows. | Fill the row or drop it. Fix "five dyes" to match. |
| 497 | (16,16) sweep | Minor | Drop it: (14,14) already shows saturation. |
| 552 | THC rank 22/44 gap error against exact CASCI | The cost model assumes THC preserves gaps to well within 0.07 eV. Untested. | Recommended. This is the one QC-side number a quantum referee will ask for. |
| 653 | "Option B": QPE + NEVPT2 with quantum RDMs | Not needed for this paper | Delete the TODO and keep the sentence "the measurement cost they add has to fit the same budget" as future work. |
| 708 | SI contents | Submission blocked | See Section 3. |
| 721 | Zenodo DOI | JCTC data policy | Required. |

Also, before submission:

- **Reference format.** Cite Zhou et al. as the PRA version (the related-work notes say it was re-checked
  against the published version; confirm that the 177 qubits, 7.96×10⁷ Toffolis and 483 shots are
  unchanged in the PRA text).
- **arXiv entries.** Eklund2026, Zak2026 and Genin2025 are arXiv entries; check them for journal versions.
- **Build comment.** Remove the comment in `main.tex` l. 4 about `main_jpca_draft.tex` (harmless, but
  clean the source before uploading).
- **Running heads.** `si.tex` repeats the full 27-word title in `\title`. Fine, but the short title should
  read "Supporting Information: Accuracy budget and quantum cost…".

---

## 5. Figures, tables and SI cross-references

**Numbering and citation**

- All 6 tables and 3 figures in the main text, and all 3 tables and 3 figures in the SI, are cited.
- No reference is undefined.
- The order of first citation is wrong for Figures 1/2 and Tables 1/2 (R6).
- The main text cites the SI correctly: Figure S1 = absorption parity, Figure S2 = active-space error,
  Figure S3 = AE demo; Tables S1–S3 = canonical CASCI/NEVPT2, SA-NEVPT2, sweep. These all match the SI
  build.
- `si.tex` l. 49 hard-codes "Table 3 of the main text". It is correct now, but it will break as soon as
  C9 or R6 is applied. Use the `xr` package (`\externaldocument{main}`) and `\ref{tab:window}`, or
  re-check this by hand after every change.

**Captions (self-containedness)**

- *Table 1.* Fine. Add "error = calc − exp".
- *Table 2.* Define ΔE_ST as the S1–T1 gap (vertical calc against experimental). Say where the
  experimental ΔE_ST comes from (phosphorescence or estimate). Say that the BDP-0 and I-BDP-0 rows are
  omitted and why. k_ISC calc for aza-BDP-Ph₄ is 1.5×10⁵ here but 1.4×10⁵ in `decision_cost.md`; make
  them consistent.
- *Table 3.* Good. Add the basis (def2-SVP) and the geometry source.
- *Table 4.* Add "gas phase", "SC-NEVPT2", "experimental values from ref Ly2021 / Killoran", and that the
  T2–S1 values are CASCI or NEVPT2 triplets on singlet-optimized orbitals.
- *Table 5.* Fix the column mismatch (R4). Say that V1 is computed with the same effective SOC operator.
  Say why the ΔE column is so far from experiment ("CASCI in CAS(4,4)/(6,6)").
- *Table 6.* Define k_th = k_r + k_ic in the caption. State the orbitals used for α_SO (R5). State
  "8 of 15 dyes; all in Table S#". State the overlap assumption γ² = 0.5 and that MLAE uses the
  empirical constant 45.
- *Figure 1.* Give the three λ_s values in the caption, and state S = 0.8 and ħω = 0.17 eV.
- *Figure 2.* Mark azaBDP\_Br2 (0.016 eV) on the plot or say in the caption that it lies outside the
  band.
- *Figure 3.* Good. Add "THC (12e,11o), 158 qubits, (1.4–2.1)×10⁸ Toffolis per call".
- *SI Figure S2.* Name the dyes and the sign convention (calc − exp).
- *SI Table S1.* The TDDFT+ΔSCF column duplicates Table 1; acceptable.

---

## 6. Headline-number consistency (abstract / results / discussion / conclusions)

Each number was checked against the tables and against `paper/quantum/results/decision_cost.md`.

| Quantity | Abstract | Results | Discussion | Conclusions | Status |
|---|---|---|---|---|---|
| Factor-2 gap tolerance | 0.06–0.08 eV | 0.062–0.083 eV (l. 372) | "about 0.07 eV" (l. 641) | "about 0.07 eV (0.03–0.13 with λ_s)" | Consistent. Use "0.06–0.08 eV" plus "≈0.07 eV" once each, and keep the λ_s range in one place. |
| Per-state accuracy | ~2 mHa | 2 mHa (l. 376) | – | – | OK (0.08/√2 eV = 2.1 mHa) |
| Top-candidate tolerance | 0.016 eV | 0.016 eV | – | 0.016 eV | OK |
| Pipeline absorption error | – | −0.055 ± 0.02 eV compact; −0.21 to −0.23 distyryl | – | – | OK against Table 1 |
| Canonical window π weight | 2.5–6.7 of 11 | 2.5–6.7 (Table 3) | – | – | OK. The abstract says "holds 2.5–6.7 of the 11 orbitals", but these are summed weights, not orbitals; say "π character equivalent to 2.5–6.7 of 11 orbitals". |
| π-space NEVPT2 S1–T1 error | 0.33–0.37 eV, "five times" | 0.33–0.37 (Table 4: 0.33, 0.35, 0.37) | "five times" | "about five times" | OK, subject to R1 |
| π-space NEVPT2 S1 error | – | 0.54–0.59 eV | – | – | OK (0.54, 0.55, 0.54, 0.59) |
| Bare CASCI S1 error | – | 1.4–2.0 eV canonical/RHF (l. 487); 1.2–2.3 eV AVAS without orbital optimization (l. 461); **0.40–0.56 eV AVAS on SA-CASSCF orbitals (Table 4, not stated)** | "1–2 eV" | "20–30 times the budget" | **Inconsistent / misleading** (R2) |
| Filter verification | "verified against exact results" | 0.03 % (2 % smallest rate) | – | "verified to 0.03 %" | Add "(2 % at the dynamic-range floor)" in the Conclusions, or say "to better than 0.1 % for all but the smallest rate". |
| Toffolis per call | ~2×10⁸ | 1.4–2.1×10⁸ | – | ~2×10⁸ | OK; see R15 |
| Logical qubits | – | 131–374; 158 used | – | 131–374 | See R16 |
| Calls per decision | 10³–10⁶ | 3×10⁴–3×10⁶ (D1); D2 down to 4×10² (not stated) | – | 10³–10⁶ | Lower end slightly off (R15) |
| Light-atom negatives (D2) | ~10¹¹ | 5×10¹⁰–5×10¹¹ | – | "between 10¹¹ and 10¹⁴" | Abstract and Conclusions round 5×10¹⁰ up to 10¹¹; acceptable as "about", but "between 10¹¹" in the Conclusions understates the minimum. Use "5×10¹⁰–4×10¹⁴". |
| Near threshold (D2) | up to 10¹⁴ | 2–4×10¹⁴ | – | 10¹⁴ | OK |
| D1 MLAE range | – | 5×10¹²–6×10¹⁴ | – | – | OK against `decision_cost.md` (5.4×10¹², 6.1×10¹⁴) |
| Sampling/MLAE ratio | – | 2×10²–2×10⁴ | – | – | OK (213–17,977) |
| α_SO/‖χ‖ | – | 5–160 (l. 560) | 5–160 (l. 658) | – | OK over all ISC-active dyes (I1-BDP 5.2 to azaBDP\_IPh2 156). Only 6.6–156 can be checked from Table 6. |
| α_SO range for "heavy-atom dyes" | – | 850–12,400 | – | – | azaBDP\_BrPh2 (Br) has 109 cm⁻¹. Say "for the iodinated and dibrominated dyes" or give the full range. |
| p range | – | 6×10⁻⁹–4×10⁻⁵ | – | – | OK |
| Versus Zhou 3.8×10¹⁰ | – | D1 10²–10⁴× | – | – | OK (142–16,000×) |
| Campaign factor N = 100 | – | 2.2 | – | – | OK (ln 4000 / ln 40 = 2.25) |
| k_ISC I2-BDP / I3-BDP | – | 6.9/7.1×10⁹ against 3.5/3.9×10⁹ | – | – | OK; "within factor 2" holds (1.97, 1.82) |
| I2-BDP-0 | – | 1.7×10⁷, 26× below; Φ_T 0.04 against 0.76 | – | – | OK (4.5×10⁸ / 1.7×10⁷ = 26) |
| Limitations: "five dyes so far" (π space) | – | Table 4 has 5 CASCI rows but only 4 NEVPT2 rows | – | – | Clarify |

---

## 7. Minor issues

- **l. 538.** The section title "Fault-tolerant resources, like for like" is colloquial. Use "Fault-tolerant
  resources compared with the gap-free proxy".
- **l. 683.** "the screening decision needs these gaps": "these" has no antecedent in the Conclusions.
  Write "singlet–triplet gaps".
- **l. 99.** "routinely computed with golden-rule approaches" is fine. The paragraph l. 78–101 is one long
  citation run. Split it after l. 94.
- **l. 164–165.** λ_s, ħω and S are "set a priori". Give the source or rationale in one clause (typical
  BODIPY C=C/C=N stretch, etc.), since the budget depends on λ_s.
- **l. 166–167.** "This is an approximation for the aza-BODIPYs, whose non-radiative decay is faster."
  Give a number and say in which direction this biases Φ_T, since three screening candidates are
  aza-BODIPYs and one of them is the "top candidate".
- **l. 179.** Absorption window "650–850 nm". Zhou et al. use 700–850 nm (`xanadu_comparison.md`). State
  the difference and why.
- **l. 204, l. 239.** File paths in the Methods (`psscreen/isc_filtered.py`, `paper/quantum/decision_cost.py`)
  belong in Data Availability only. They already appear there.
- **l. 218.** "and K′ is known exactly from the Hamiltonian" is useful, but it is buried at the end of a
  long sentence. Make it its own sentence: this is why D2 can sometimes be decided without any call.
- **l. 226.** The degree formula carries a −0.11 offset term inside parentheses multiplied by λ′/Δ;
  `decision_cost.md` prints "+ -0.106". Check the sign and print it as "− 0.11".
- **l. 319–321.** Good robustness check. Also give the factor relative to experiment at the extremes
  (2.3×10⁹ is 0.7× experiment, 2.8×10¹⁰ is 8×).
- **l. 395.** "for gap errors of 0.01–0.02 eV" are smaller than the top-candidate tolerance of 0.016 eV,
  yet the paragraph title says the top candidate is "the most fragile". Clarify: most fragile in k_ISC;
  the near-threshold dyes are the most fragile in class.
- **l. 506.** "to all printed digits": see R4.
- **l. 591.** `azaBDP\_IPh2` and `Br2-BDP` are labels from the code (R7).
- **l. 595–597.** "D2 for light-atom negatives is comparable to it". Comparable to 3.8×10¹⁰: ours are
  5×10¹⁰–5×10¹¹, so "1–10 times" is more precise.
- **l. 669–678 (Limitations).** Good and honest. Add "no solvent in the active-space references either"
  and "no basis-set convergence for NEVPT2" (R1).
- **Units.** Use "cm⁻¹" consistently. In Table 5, V1 "<10⁻³" and "1×10⁻³" have different precision.
- **Acknowledgement (l. 702).** Thanking "the Xanadu team, authors of ref" is gracious. Check that they
  agree to be named, or say "the authors of ref X".
- **Hyphenation.** Mixed use of "Spin–orbit", "spin-orbital" (l. 211) and "spin–orbit-active" (l. 697).
  Use "spin-orbital" only for the orbital basis.
