"""El hueco del backend cuántico: CASCI exacto vs un 'solver externo' (aquí FCI de PySCF
invocado a través de la interfaz de callable) deben coincidir; se exporta FCIDUMP."""
import numpy as np, sys, os, tempfile
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from psscreen.backends import ActiveSpaceBackend

ATOM = "C 0 0 0\nO 0 0 1.205\nH 0 0.94 -0.587\nH 0 -0.94 -0.587"

def external_solver(h1, h2, ncas, nelecas, nroots, spin):
    """Ejemplo de solver externo: aquí iría QPE/QSP. Firma: devuelve energías activas (Ha)."""
    from pyscf import fci
    s = fci.direct_spin1.FCI(); s.nroots = nroots
    fci.addons.fix_spin_(s, ss=spin / 2 * (spin / 2 + 1))
    e, _ = s.kernel(h1, h2, ncas, nelecas)
    return np.atleast_1d(e)

def test_quantum_slot():
    d = tempfile.mkdtemp()
    a = ActiveSpaceBackend(ncas=4, nelecas=4, basis="6-31g", nroots=2,
                           fcidump_path=os.path.join(d, "FCIDUMP")).run(ATOM)
    b = ActiveSpaceBackend(ncas=4, nelecas=4, basis="6-31g", nroots=2, solver=external_solver).run(ATOM)
    print("CASCI S:", a.e_singlets, "T:", a.e_triplets)
    assert os.path.exists(os.path.join(d, "FCIDUMP"))
    assert np.allclose(a.e_singlets, b.e_singlets, atol=1e-6) and np.allclose(a.e_triplets, b.e_triplets, atol=1e-6)

if __name__ == "__main__":
    test_quantum_slot(); print("OK")
