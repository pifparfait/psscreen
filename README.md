# psscreen: photosensitizer screening (absorption + ISC)

A framework for ranking candidate photosensitizers for photodynamic therapy (PDT) on two criteria:
absorption in the therapeutic window and S1→Tn intersystem crossing (ISC). It starts from the
approach of arXiv:2512.15889 (Xanadu), fixes its weak points, and provides a classical backend
validated against experiment plus a PennyLane quantum backend that implements the paper's two
algorithms on small active spaces (simulator).

## Pipeline

```
SMILES ──RDKit/MMFF──► GFN2-xTB (tblite) ──► B3LYP/def2-SVP (PySCF, density fitting)
        │
        ├─ Singlets: TDDFT (RPA), own Davidson solver, full space         psscreen/tdwin.py
        ├─ Triplets: explicit TDA in an orbital window                    psscreen/tdwin.py
        ├─ ΔSCF/ROKS correction of S1 and T1 (the "cyanine problem")      psscreen/deltascf.py
        ├─ SOC: effective 1e operator calibrated on atomic spectra        psscreen/soc.py
        │        <S_m|H_SO|T_n,M> from TDA amplitudes, summed over M
        ├─ Rates: k_r (Einstein/Strickler-Berg), k_ISC golden rule with   psscreen/rates.py
        │        Marcus-Levich-Jortner Franck-Condon factors (S-T gap
        │        included), Φ_T, Φ_F
        └─ PDT score: window × brightness × Φ_T × (T1 > O2 1Δg)           psscreen/pipeline.py

Quantum backend (psscreen/quantum_pl.py, PennyLane):
  CAS(n,n) → 2n qubits (Jordan-Wigner)
  ├─ Threshold projection with QSVT (qml.QSVT + block encoding + PCPhase) → dipole weight in a window
  ├─ ISC proxy: |<T,M| e^{-iH_SOC t} |S1>|² with qml.TrotterProduct
  ├─ Z2 tapering (qml.taper) and FTQC resources (double factorization)
  └─ ActiveSpaceBackend (backends.py): exports FCIDUMP or accepts any external solver
```

## What it fixes relative to the paper

| Weakness of the paper | Here |
|---|---|
| No dynamic correlation (CASCI only) | Classical TDDFT + ΔSCF backend; optional CASCI/NEVPT2; the quantum result is compared against both |
| ISC proxy ≈ \|⟨T\|H_SOC\|S⟩\|², no gap or vibronic effects | k_ISC from the golden rule with Marcus-Levich-Jortner FCWD and all T_n, not only T1 |
| One-electron SOC without screening | Per-element constants calibrated to reproduce experimental atomic ζ_np (equivalent to a Koseki-type Z_eff; C 3.5, N 4.4, F 6.4 come out on their own) |
| No experimental validation | Validation against 4 molecules with published data (`validation/`) |

## Installation

### macOS (Apple Silicon), recommended: conda/mamba

```bash
# if you do not have conda: brew install miniforge
git clone git@github.com:pifparfait/psscreen.git && cd psscreen
conda env create -f environment.yml
conda activate psscreen
python tests/test_tdwin.py && python tests/test_soc.py   # quick check (~1 min)
```

Mac notes: PySCF uses OpenMP; set the thread count with `export OMP_NUM_THREADS=8` (performance
cores). If `pip install pyscf` fails on arm64, use the conda-forge build (above).

### pip only (Linux, or Mac with compilers)

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

## Usage

```bash
python screen.py candidates.csv --out runs/        # columns: name,smiles
python validation/run_validation.py               # classical (all 4 molecules, or pass names)
python validation/run_quantum.py BDP 4 6          # quantum: CAS(4,4) and CAS(6,6)
python validation/report.py                       # predicted vs experimental tables
python tests/test_tdwin.py && python tests/test_soc.py && python tests/test_quantum_pl.py && python tests/test_active_space.py
```

A real quantum backend (hardware or QPE) plugs in as `ActiveSpaceBackend(solver=f)` with
`f(h1, h2, ncas, nelecas, nroots, spin) -> energies`, or by replacing the functions in
`quantum_pl.py`, which already work on the qubit Hamiltonian.

## Validation set

| Molecule | Experiment | Reference |
|---|---|---|
| BODIPY core | λ 503 nm, Φ_F ≈ 1 | Schmitt et al., J. Fluoresc. 2009 |
| BDP (1,3,5,7-tetramethyl-8-phenyl) | λ 503 nm, Φ_F 0.45, no ISC (toluene) | Ly et al., PCCP 2021 |
| I2-BDP (2,6-diiodo) | λ 537 nm, k_ISC 3.5×10⁹ s⁻¹, Φ_T 0.84 (toluene) | Ly et al., PCCP 2021 |
| Tetraphenyl aza-BODIPY | λ ≈ 650-655 nm, Φ_F 0.34 | Killoran 2002; Gorman 2004; Loudet 2008 |

## Empirical parameters (state them when you use the framework)

- `k_ic` (S1→S0 internal conversion): 1.5e8 s⁻¹ by default (non-radiative rate measured for BDP).
- MLJ: λ_s = 0.10 eV, ħω = 0.17 eV, S = 0.8. A sensitivity analysis is included in `validation/report.py`.
- Vertical energies at the S0 geometry, gas phase (no PCM in the in-house TDDFT yet).

## Cost

A 43-atom molecule on 4 cores takes on the order of 10-45 min; 60 atoms about 1 h or more.
