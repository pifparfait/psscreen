# Paper draft

Manuscript draft for *The Journal of Physical Chemistry A* (ACS `achemso` class).

- `main.tex`: manuscript. Items in red (`\todo{...}`) need data or a decision.
- `refs.bib`: bibliography. Check every entry before submission.
- `make_figures.py`: builds `figures/*.pdf` from `validation/runs/*.json`. Run from the repository root:
  `python paper/make_figures.py`. Molecules without a result file are skipped, so the styryl dyes appear
  once their recalculated (trans) results are in `validation/runs/`.

Build: `latexmk -pdf main.tex` (needs a TeX distribution with `achemso`, e.g. TeX Live or Overleaf).

## Open items

1. Styryl dyes (BDP-0, I-BDP-0, I2-BDP-0): results with trans geometries are being recalculated; fill
   Table 1, Table 2 and the "Styryl-extended dyes" section.
2. Re-check experimental values against Ly et al. (PCCP 2021) Table 2.
3. Source of the logical-qubit difference with arXiv:2512.15889 (601 vs 180-350).
4. MLJ sensitivity table and remaining SI material.
5. Corresponding e-mail, funding, Zenodo DOI.
