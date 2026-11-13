"""Physics tests of the SOC module: El-Sayed rules in formaldehyde and rotational invariance."""
import numpy as np, sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from psscreen.backends import PySCFTDA
from psscreen.soc import soc_mo_integrals, soc_singlet_triplet, soc_ground_triplet

H2CO = [("C", [0, 0, 0]), ("O", [0, 0, 1.205]), ("H", [0, 0.94, -0.587]), ("H", [0, -0.94, -0.587])]

def run(rot=np.eye(3)):
    atom = "\n".join(f"{s} {' '.join(map(str, rot @ np.array(x)))}" for s, x in H2CO)
    be = PySCFTDA(xc="b3lyp", basis="def2-svp", nstates=3, density_fit=False)
    r = be.run(atom)
    h, _ = soc_mo_integrals(r.mol, r.mo_coeff, be.basis)
    return r, soc_singlet_triplet(r, h), soc_ground_triplet(r, h)

def test_el_sayed_and_rotation():
    r, V, V0 = run()
    print("S:", r.e_singlets.round(2), "T:", r.e_triplets.round(2))
    print("V(S1,Tn) cm-1:", V[0].round(2)); print("V(S0,Tn) cm-1:", V0.round(2))
    # S1 y T1 son n->pi*: acoplamiento directo ~0 (El-Sayed); S0-T1(n pi*) no nulo
    assert V[0, 0] < 1.0
    assert V0[0] > 10.0
    th = 0.7; R = np.array([[1, 0, 0], [0, np.cos(th), -np.sin(th)], [0, np.sin(th), np.cos(th)]]) @ \
        np.array([[np.cos(1.1), 0, np.sin(1.1)], [0, 1, 0], [-np.sin(1.1), 0, np.cos(1.1)]])
    _, V2, V02 = run(R)
    # lowest states only: Davidson may capture different high roots after rotation
    assert np.allclose(V[:2, :2], V2[:2, :2], atol=0.05) and np.allclose(V0[:2], V02[:2], atol=0.05)

if __name__ == "__main__":
    test_el_sayed_and_rotation(); print("OK")
