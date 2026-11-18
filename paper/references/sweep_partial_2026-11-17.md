# Active-space sweep, partial (from the tars logs, 2026-11-17 12:01 UTC)

Vertical energies (eV), def2-SVP, SA-CASSCF(n,n) orbitals (3 singlets), warm-started from the previous size.

| Dye | (n,n) | CASCI S1 | CASCI S1−T1 | CASCI T2−S1 | NEVPT2 S1 | NEVPT2 S1−T1 | NEVPT2 T2−S1 |
|---|---|---|---|---|---|---|---|
| BDP | (8,8) | 3.801 | 1.736 | −0.216 | 3.082 | 0.998 | +0.428 |
| BDP | (10,10) | 3.222 | 1.169 | +0.192 | 3.015 | 0.993 | +0.446 |
| BDP | (12,12) | 3.059 | 1.033 | +0.292 | 2.970 | 0.955 | +0.475 |
| I2-BDP | (8,8) | 3.161 | 1.082 | +0.157 | 2.842 | 0.816 | +0.400 |
| I2-BDP | (10,10) | 2.944 | 0.940 | +0.740 | 2.858 | 0.874 | +0.641 |
| I2-BDP | (12,12) | 2.937 | 1.070 | +0.049 | 2.836 | 0.882 | +0.522 |

Experiment: BDP S1 2.46, S1−T1 0.77; I2-BDP S1 2.31, S1−T1 0.59 (Ly et al. 2021).

Update 2026-11-18 00:18 UTC: (14,14) CASCI on SA-CASSCF(14,14) orbitals (no NEVPT2 by design):
BDP S1 3.060, S2 4.181, T1 2.028, T2 3.354, T3 4.203 (S1−T1 1.032, T2−S1 +0.294);
I2-BDP S1 2.938, S2 3.771, T1 1.869, T2 2.988, T3 4.293 (S1−T1 1.069, T2−S1 +0.050). Both within 0.003 eV of (12,12).
