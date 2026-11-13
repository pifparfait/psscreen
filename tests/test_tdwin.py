"""El TDDFT en ventana reproduce PySCF completo cuando la ventana abarca todo el espacio."""
import numpy as np, sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pyscf import gto, dft, tdscf
from psscreen.tdwin import window_td

def test_matches_pyscf():
    mol = gto.M(atom="C 0 0 0; O 0 0 1.205; H 0 0.94 -0.587; H 0 -0.94 -0.587", basis="6-31g*", verbose=0)
    mf = dft.RKS(mol, xc="b3lyp").density_fit().run()
    r = window_td(mf, ewin_ev=1e4, nstates=4, max_occ=999, max_vir=999)
    ref = tdscf.TDDFT(mf); ref.nstates = 4; ref.kernel()
    refT = tdscf.TDA(mf); refT.singlet = False; refT.nstates = 4; refT.kernel()
    print(r["e_S"] * 27.2114, ref.e * 27.2114); print(r["f"], ref.oscillator_strength())
    print(r["e_T"] * 27.2114, refT.e * 27.2114)
    assert np.allclose(r["e_S"], ref.e, atol=2e-4)
    assert np.allclose(r["f"], ref.oscillator_strength(), atol=2e-3)
    assert np.allclose(r["e_T"], refT.e, atol=2e-4)



def test_rpa_davidson_matches_pyscf():
    from psscreen.tdwin import rpa_singlets_davidson
    mol = gto.M(atom="C 0 0 0; O 0 0 1.205; H 0 0.94 -0.587; H 0 -0.94 -0.587", basis="6-31g*", verbose=0)
    mf = dft.RKS(mol, xc="b3lyp").density_fit().run()
    r = rpa_singlets_davidson(mf, nstates=3, ewin_ev=1e4, grid_level=3, tol=1e-6)
    # referencia: diagonalización exacta del espacio completo (el Davidson de PySCF
    # a veces se salta una raíz casi degenerada, como S3 aquí)
    ref = window_td(mf, ewin_ev=1e4, nstates=3, max_occ=999, max_vir=999)
    print(r["e_S"] * 27.2114, ref["e_S"] * 27.2114, r["f"], ref["f"])
    assert np.allclose(r["e_S"], ref["e_S"], atol=1e-5)
    assert np.allclose(r["f"], ref["f"], atol=1e-4)

if __name__ == "__main__":
    test_matches_pyscf(); test_rpa_davidson_matches_pyscf(); print("OK")
