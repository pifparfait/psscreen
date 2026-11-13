"""SMILES -> geometría 3D (RDKit ETKDG/MMFF) -> optimización GFN2-xTB (tblite + ASE)."""
from __future__ import annotations
import numpy as np


def smiles_to_xyz(smiles: str, seed: int = 7) -> tuple[list[str], np.ndarray]:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if AllChem.EmbedMolecule(mol, randomSeed=seed) != 0:
        raise RuntimeError(f"Embedding falló para {smiles}")
    if AllChem.MMFFHasAllMoleculeParams(mol):
        AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
    else:
        AllChem.UFFOptimizeMolecule(mol, maxIters=2000)
    sym = [a.GetSymbol() for a in mol.GetAtoms()]
    xyz = mol.GetConformer().GetPositions()
    return sym, np.asarray(xyz)


def xtb_optimize(sym, xyz, charge: int = 0, fmax: float = 0.02, steps: int = 500):
    """Optimiza con GFN2-xTB. Devuelve (sym, xyz_Å, energía_eV)."""
    from ase import Atoms
    from ase.optimize import BFGS
    from tblite.ase import TBLite
    atoms = Atoms(symbols=sym, positions=xyz)
    atoms.calc = TBLite(method="GFN2-xTB", charge=charge, verbosity=0)
    BFGS(atoms, logfile=None).run(fmax=fmax, steps=steps)
    return sym, atoms.get_positions(), atoms.get_potential_energy()


def to_pyscf_atom(sym, xyz) -> str:
    return "\n".join(f"{s} {x:.6f} {y:.6f} {z:.6f}" for s, (x, y, z) in zip(sym, xyz))


def write_xyz(path, sym, xyz, comment=""):
    with open(path, "w") as f:
        f.write(f"{len(sym)}\n{comment}\n")
        for s, (x, y, z) in zip(sym, xyz):
            f.write(f"{s} {x:.6f} {y:.6f} {z:.6f}\n")


def read_xyz(path):
    lines = open(path).read().splitlines()
    n = int(lines[0])
    sym, xyz = [], []
    for l in lines[2:2 + n]:
        p = l.split()
        sym.append(p[0]); xyz.append([float(v) for v in p[1:4]])
    return sym, np.array(xyz)
