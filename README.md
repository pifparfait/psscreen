# psscreen: cribado de fotosensibilizadores (absorción + ISC)

Framework para clasificar candidatos a fotosensibilizador para terapia fotodinámica (PDT)
con dos criterios: absorción en la ventana terapéutica e intersystem crossing (ISC) S1→Tn.
Parte del planteamiento de arXiv:2512.15889 (Xanadu) y corrige sus puntos débiles, con un
backend clásico validado contra experimento y un backend cuántico en PennyLane que implementa
los dos algoritmos del paper sobre espacios activos pequeños (simulador).

## Pipeline

```
SMILES ──RDKit/MMFF──► GFN2-xTB (tblite) ──► B3LYP/def2-SVP (PySCF, DF)
        │
        ├─ Singletes: TDDFT (RPA) con Davidson propio, espacio completo   psscreen/tdwin.py
        ├─ Tripletes: TDA explícito en ventana de orbitales               psscreen/tdwin.py
        ├─ Corrección ΔSCF/ROKS de S1 y T1 (problema de las cianinas)     psscreen/deltascf.py
        ├─ SOC: operador efectivo 1e calibrado con espectros atómicos     psscreen/soc.py
        │        ⟨S_m|H_SO|T_n,M⟩ con amplitudes TDA, suma sobre M
        ├─ Tasas: k_r (Einstein/Strickler-Berg), k_ISC regla de oro con    psscreen/rates.py
        │        Franck-Condon Marcus-Levich-Jortner (gap S-T incluido), Φ_T, Φ_F
        └─ Puntuación PDT: ventana × brillo × Φ_T × (T1 > 1Δg O2)          psscreen/pipeline.py

Backend cuántico (psscreen/quantum_pl.py, PennyLane):
  CAS(n,n) → 2n qubits (Jordan-Wigner)
  ├─ Proyección por umbral con QSVT (qml.QSVT + BlockEncode + PCPhase) → peso dipolar en ventana
  ├─ Proxy ISC: |⟨T,M| e^{-iH_SOC t} |S1⟩|² con qml.TrotterProduct
  ├─ Tapering Z2 (qml.taper) y recursos FTQC (qml.resource.DoubleFactorization)
  └─ ActiveSpaceBackend (backends.py): exporta FCIDUMP o acepta cualquier solver externo
```

## Qué corrige respecto al paper

| Debilidad del paper | Aquí |
|---|---|
| Sin correlación dinámica (CASCI) | Backend clásico TDDFT + ΔSCF; CASCI/NEVPT2 opcional; el resultado cuántico se compara con ambos |
| Proxy ISC ≈ \|⟨T\|H_SOC\|S⟩\|², sin gap ni vibraciones | k_ISC por regla de oro con FCWD de Marcus-Levich-Jortner y todos los T_n, no solo T1 |
| SOC de un electrón sin apantallamiento | Constantes por elemento calibradas para reproducir ζ_np atómicos experimentales (equivale a Z_eff de Koseki; C 3.5, N 4.4, F 6.4 salen solos) |
| Sin validación experimental | Validación contra 4 moléculas con datos publicados (validation/) |

## Instalación

### macOS (Apple Silicon) — recomendado con conda/mamba

```bash
# si no tienes conda: brew install miniforge
git clone git@github.com:pifparfait/psscreen.git && cd psscreen
conda env create -f environment.yml
conda activate psscreen
python tests/test_tdwin.py && python tests/test_soc.py   # comprobación rápida (~1 min)
```

Notas Mac: PySCF usa OpenMP; fija los hilos con `export OMP_NUM_THREADS=8` (núcleos de
rendimiento). Si `pip install pyscf` falla en arm64, usa la versión de conda-forge (arriba).

### Solo pip (Linux o Mac con compiladores)

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

## Uso

```bash
python screen.py candidatos.csv --out runs/        # columnas: name,smiles
python validation/run_validation.py               # clásico (las 4 moléculas; o pasa nombres)
python validation/run_quantum.py BDP 4 6          # cuántico, CAS(4,4) y CAS(6,6)
python validation/report.py                       # tabla predicho vs experimento
python tests/test_tdwin.py && python tests/test_soc.py && python tests/test_quantum_pl.py && python tests/test_active_space.py
```

Un backend cuántico real (hardware o QPE) se enchufa como `ActiveSpaceBackend(solver=f)` con
`f(h1, h2, ncas, nelecas, nroots, spin) -> energías`, o reemplazando las funciones de
`quantum_pl.py`, que ya trabajan con el Hamiltoniano de qubits.

## Parámetros empíricos (decláralos al usarlo)

- `k_ic` (conversión interna S1→S0): 1.5e8 s⁻¹ por defecto (k_nr medido del BDP).
- MLJ: λ_s = 0.10 eV, ħω = 0.17 eV, S = 0.8. Se incluye análisis de sensibilidad.
- Energías verticales en la geometría S0, fase gas (sin PCM todavía en el TDDFT propio).

## Coste

Molécula de 43 átomos con 4 núcleos: del orden de 10-30 min. 60 átomos: ~1 h.
