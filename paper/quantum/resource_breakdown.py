"""Break down the PennyLane double-factorization QPE estimate for CAS(12,11) into its qubit registers.

Explains why our logical-qubit count (about 600) is larger than the 177 qubits that arXiv:2512.15889
reports for the threshold projection of the same active-space size: most of our qubits are the QROM output
register for the Givens-rotation angles (n * beta * kr / 2, Lee et al., PRX Quantum 2, 030305, Eq. C40).
Usage (from the repository root): python paper/quantum/resource_breakdown.py [NAME]
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from pennylane.estimator.qpe_resources.second_quantization import DoubleFactorization as DF  # noqa: E402
from pyscf import dft, gto  # noqa: E402

from psscreen.backends import _ecp  # noqa: E402
from psscreen.geometry import read_xyz, to_pyscf_atom  # noqa: E402
from psscreen.quantum_pl import active_space  # noqa: E402

name = sys.argv[1] if len(sys.argv) > 1 else "BODIPY_core"
sym, xyz = read_xyz(os.path.join(ROOT, "validation", "runs", f"{name}.xyz"))
atom = to_pyscf_atom(sym, xyz)
mol = gto.M(atom=atom, basis="def2-svp", ecp=_ecp(atom, "def2-svp"), verbose=0)
mf = dft.RKS(mol, xc="b3lyp").density_fit()
mf.grids.level = 2
mf.kernel()
asp = active_space(mf, 11, 12, "def2-svp")
n = 2 * asp["ncas"]
print("| settings | logical qubits | system qubits | rotation-angle register | Toffolis | lambda (Ha) |")
print("|---|---|---|---|---|---|")
for kw in ({}, {"beta": 16}, {"beta": 12}, {"error": 0.0043}):
    df = DF(asp["h1"], asp["h2"].transpose(0, 2, 3, 1), **kw)
    beta = kw.get("beta", 20)
    rr = df.rank_r * df.rank_m
    kr = DF._qrom_cost((rr, n / 2, rr, n * beta, 0))[1]
    label = ", ".join(f"{k}={v}" for k, v in kw.items()) or "default (error 1.6 mHa, beta=20)"
    print(f"| {label} | {df.qubits} | {n} | {int(kr * n * beta / 2)} (kr={kr}) | {df.gates:.2e} | {df.lamb:.1f} |")
