"""Total logical qubits of one decision circuit (M9 of the mock review), with PennyLane's estimator.

Walk: THC qubitization of the (12e,11o) BODIPY-core Hamiltonian, rank 33, angles loaded in five batches, the
precisions of thc_pennylane_resources.py. Filter: QSVT of that walk at the degree of decision_cost.json. Canonical
amplitude estimation: phase estimation on the filter with enough phase bits for the smallest relative error needed
(D1: r = 0.2 at the smallest p of the ISC-active dyes). The one-body H_SO block encoding (a fixed basis rotation of
the 22 spin orbitals, a 5-qubit index register prepared with |mu_k| weights, one sign qubit) and the S1 eigenstate
filter (one QSVT signal qubit) are added by hand; they act before or after the walk and could reuse its ancillas,
so the sum is an upper bound. Usage (repository root): python paper/quantum/logical_qubits.py
Writes paper/quantum/results/logical_qubits.json.
"""
import json
import math
import os

import numpy as np
import pennylane.estimator as qre

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
N, M, LAM, EPS = 11, 33, 15.9, 0.16e-3
ALEPH = int(np.ceil(2.5 + np.log2(LAM / EPS)))
BETH = int(np.ceil(5.625 + np.log2(2 * LAM * N / EPS)))
dc = json.load(open(os.path.join(ROOT, "paper", "quantum", "results", "decision_cost.json")))


def wires(op):
    return int(str(qre.estimate(op)).splitlines()[1].split(":")[1].strip())


def walk():
    h = qre.THCHamiltonian(num_orbitals=N, tensor_rank=M, one_norm=LAM)
    return qre.QubitizeTHC(h, prep_op=qre.PrepTHC(h, coeff_precision=ALEPH),
                           select_op=qre.SelectTHC(h, num_batches=5, rotation_precision=BETH))


active = [r for r in dc["rows"] if r.get("k_isc", 0) >= 1e6]
deg = int(max(r["D1"]["degree"] for r in active))
pmin = min(r["p"] for r in active)
t_phase = math.ceil(math.log2(math.pi / (dc["model"]["r_D1"] * math.sqrt(pmin))))   # phase bits, theta ~ sqrt(p)
w_walk = wires(walk())
w_filter = wires(qre.QSVT(walk(), encoding_dims=2 ** (2 * N), poly_deg=deg))
so_index, so_sign, s1_signal, real_part = 5, 1, 1, 1
total_mlae = w_filter + real_part + so_index + so_sign + s1_signal
total_canonical = total_mlae + t_phase
out = {"walk_wires": w_walk, "filter_degree": deg, "filter_wires": w_filter, "real_part_lcu": real_part,
       "H_SO_index": so_index, "H_SO_sign": so_sign, "S1_filter_signal": s1_signal, "p_min": pmin,
       "canonical_phase_bits": t_phase, "total_MLAE": total_mlae, "total_canonical": total_canonical}
json.dump(out, open(os.path.join(ROOT, "paper", "quantum", "results", "logical_qubits.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
