"""Styryl dyes must be built as E isomers: SMILES without E/Z labels gave Z styryl groups in an earlier run."""
import os, sys
import numpy as np
import pytest
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(here), "validation")); sys.path.insert(0, os.path.dirname(here))
Chem = pytest.importorskip("rdkit.Chem")
from molecules import SMILES
from psscreen.geometry import smiles_to_xyz

PATT = Chem.MolFromSmarts("[#6;R]-[CH]=[CH]-c")


def _dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / np.linalg.norm(b1)
    v, w = b0 - b1n * (b0 @ b1n), b2 - b1n * (b2 @ b1n)
    return np.degrees(np.arctan2(np.cross(b1n, v) @ w, v @ w))


@pytest.mark.parametrize("name", [n for n, s in SMILES.items() if "/C=C/" in s])
def test_styryl_is_trans(name):
    mol = Chem.AddHs(Chem.MolFromSmiles(SMILES[name]))
    _, xyz = smiles_to_xyz(SMILES[name])
    matches = mol.GetSubstructMatches(PATT)
    assert len(matches) == 2, "expected two styryl groups"
    for m in matches:
        assert abs(_dihedral(*[xyz[i] for i in m])) > 150, f"{name}: cis styryl"
