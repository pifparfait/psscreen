# SA-CASSCF(12,11)/SC-NEVPT2 results, partial (copied from the tars logs on 2026-11-16)

Vertical energies in eV relative to S0, def2-SVP, gas phase. The JSON files will arrive on branch
results-tars-paper when all nine jobs finish; this file records the three finished dyes so the paper can use them.

| Dye | S1 | S2 | T1 | T2 | T3 | S1 exp | S1 − T1 (exp) | T2 − S1 |
|---|---|---|---|---|---|---|---|---|
| BODIPY_core | 3.117 | 4.587 | 1.958 | 4.031 | 4.633 | 2.46 | 1.159 (—) | +0.914 |
| BDP | 3.048 | 4.030 | 2.059 | 3.472 | 3.920 | 2.46 | 0.989 (0.77) | +0.424 |
| I-BDP | 2.993 | 3.893 | 1.962 | 3.417 | 4.654 | 2.46 | 1.031 (0.76) | +0.424 |

TDDFT/ΔSCF pipeline for comparison: S1 2.45/2.40/2.39 eV, S1 − T1 0.71/0.72/0.72 eV, T2 − S1 +0.44/+0.48/+0.48 eV.

Update 2026-11-16 23:38 UTC (tars log): I3-BDP S1 2.829, S2 3.856, T1 2.030, T2 3.280, T3 3.621 eV;
S1 exp 2.30; S1 − T1 0.799 (exp 0.58); T2 − S1 +0.451 (pipeline +0.31).
