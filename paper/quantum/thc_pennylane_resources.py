"""Qubits and Toffolis per walk step for THC qubitization of CAS(12,11) with PennyLane's estimator,
to compare with arXiv:2512.15889 (THC block encoding, 177 qubits for N = 11).

The qubit count is a space-time trade-off: loading the Givens-rotation angles in `num_batches` batches
shrinks the temporary angle register at a small Toffoli cost. Precisions follow the preprint:
coefficient bits = ceil(2.5 + log2(lambda/eps)), rotation bits = ceil(5.625 + log2(2 lambda N/eps)), eps = 0.16 mHa.
lambda = 15.9 Ha is the double-factorization one-norm of the BODIPY-core CAS(12,11) Hamiltonian
(paper/quantum/resource_breakdown.py); the THC one-norm is of the same size.
Usage: python paper/quantum/thc_pennylane_resources.py
"""
import numpy as np
import pennylane.estimator as qre

N, LAM, EPS = 11, 15.9, 0.16e-3
ALEPH = int(np.ceil(2.5 + np.log2(LAM / EPS)))
BETH = int(np.ceil(5.625 + np.log2(2 * LAM * N / EPS)))

print(f"coefficient bits {ALEPH}, rotation bits {BETH}")
print("| THC rank | angle batches | QROM swap depth | total wires | Toffolis per walk step |")
print("|---|---|---|---|---|")
for M in (2 * N, 4 * N):
    for nb in (1, 2, 5, 10):
        for ssd in (None, 1):
            h = qre.THCHamiltonian(num_orbitals=N, tensor_rank=M, one_norm=LAM)
            prep = qre.PrepTHC(h, coeff_precision=ALEPH, select_swap_depth=ssd)
            sel = qre.SelectTHC(h, num_batches=nb, rotation_precision=BETH, select_swap_depth=ssd)
            res = qre.estimate(qre.QubitizeTHC(h, prep_op=prep, select_op=sel))
            lines = str(res).splitlines()
            wires = lines[1].split(":")[1].strip()
            tof = [x.split(":")[1].strip(" ,") for x in lines if "Toffoli" in x][0]
            print(f"| {M} | {nb} | {ssd or 'optimal'} | {wires} | {tof} |")
