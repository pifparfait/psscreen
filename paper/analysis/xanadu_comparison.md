# What arXiv:2512.15889 does and does not claim, and how our numbers compare (item 6 of the paper plan)

Source: Zhou, Casares, Dhawan, Loaiza, Jahangiri, Lang, Arrazola, Fomichev, "Quantum Algorithms for
Photoreactivity in Cancer-Targeted Photosensitizers", arXiv:2512.15889v1 (17 Dec 2025). Read from the arXiv
HTML version; re-check every point against the PDF and any later version before citing.

## What the preprint does

- **Absorption.** Threshold projection with two Heaviside-type QSVT polynomials (pyqsp phases) on a qubitized
  THC block encoding (BLISS-THC cost model), with the window set to 700-850 nm. The observable is the
  normalized dipole weight in the window. The polynomial degree is fitted empirically as
  d = 4.7571 λ'/Δ + 321.2 at ε_H ≈ 0.01, where λ' = λ + |E_th|.
- **ISC.** The evolution proxy |⟨T1,M| e^{-i t H_SOC} |S1⟩|², with a one-body H_SOC that is fast-forwarded by a
  basis change, so no block encoding is needed. It covers **S1→T1 only**, has **no energy gap and no
  Franck-Condon factor**, and the authors state that it "is not an absolute ISC rate". It is meant to be
  calibrated against a reference molecule. A vibronic algorithm is mentioned for the weak-SOC regime and
  estimated at about 10¹⁰ Toffolis, but it is not developed.
- **Systems.** Four BODIPY-type molecules (Orig-BDP, a triazole aza-BODIPY, Br-aza-BDP, Pt-BDP), in active
  spaces of 11-45 orbitals. The minimal one is the CAS(12,11) of Momeni and Brown (JCTC 2015). The preprint
  says CASPT2 with that active space agrees closely with experiment, and it includes a frozen PCM in the
  integrals.
- **Resources (Table I, per shot).** Orig-BDP N = 11: 177 qubits and 2.72×10⁷ Toffolis for the threshold
  projection, and 176 qubits and 7.96×10⁷ Toffolis for the ISC proxy. Up to N = 45: about 350 qubits and
  2.5-5.8×10⁹ Toffolis.
- **Not included.** The preprint has no computed spectra, no ISC values and no comparison with experiment
  or with a classical method for the quantum observables.

There are small internal inconsistencies in the text, for example "fewer than 350 qubits" while the table
lists 353. These are not worth raising in our paper.

## How to state our comparison fairly

1. **Do not attack a claim the preprint does not make.** The preprint does not claim that CASCI energies are
   accurate. It cites CASPT2 for accuracy and computes CASCI-level spectra. Our point is narrower: the quantum
   observables are evaluated on the bare active-space Hamiltonian, and at CAS(12,11) that Hamiltonian misses S1
   by 1.4-2.0 eV. A dynamic-correlation step is therefore required, and the preprint does not include one in
   its cost or workflow.
2. **Our CAS(4,4)/(6,6) simulations test the algorithms, not the preprint's accuracy.** Say this explicitly.
3. **The ISC proxy is T1-only and gap-free by design.** Our contributions are the accuracy budget, which shows
   when that matters (an upper triplet near resonance, as in azaBDP_Br2 and probably I2-BDP-0), and the
   gap-resolved filter (`psscreen/isc_filtered.py`). The filter gives an absolute golden-rule rate over all
   triplets, given a vibronic line-shape model.

## Qubit-count discrepancy: explained

Our "601 logical qubits" and the preprint's "177 qubits" do not describe the same algorithm:

| | ours | arXiv:2512.15889 |
|---|---|---|
| algorithm | QPE of the energy (qubitization) | threshold projection (QSVT) of the dipole weight |
| block encoding | double factorization (PennyLane `DoubleFactorization`, Lee et al. 2021) | BLISS-THC |
| qubits for CAS(12,11) | 601 = 22 system + 440 rotation-angle QROM register + other ancillas | 177 = 22 system + max(n_SoS, n_aux + 3) |

In the double-factorization model, 440 of the 601 qubits are the QROM output register for the Givens
rotation angles: n·β·k_r/2 with n = 22 spin-orbitals, β = 20 bits and k_r = 2, a space-time trade-off from
Eq. (C40) of Lee et al. Lowering β to 12 bits gives 417 qubits, and the Toffoli count changes from 5.7×10⁷ to
4.3×10⁷ (BODIPY core, `paper/quantum/resource_breakdown.py`). The Toffoli counts are of the same order in both
works (2-6×10⁷ per run or shot for N = 11). In the paper we should state the two algorithms separately and
not present the qubit difference as a discrepancy.

Values for the BODIPY core, CAS(12,11), from `resource_breakdown.py`:

| settings | logical qubits | rotation-angle register | Toffolis | λ (Ha) |
|---|---|---|---|---|
| default (1.6 mHa, β = 20) | 601 | 440 | 5.7×10⁷ | 15.9 |
| β = 16 | 509 | 352 | 5.0×10⁷ | 15.9 |
| β = 12 | 417 | 264 | 4.3×10⁷ | 15.9 |
| error 4.3 mHa | 599 | 440 | 2.1×10⁷ | 15.9 |
