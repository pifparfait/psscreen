"""Conformer search for flexible dyes: RDKit ETKDG -> GFN2-xTB optimization -> energy/RMSD filtering.

Writes the retained conformers as <out>/<NAME>_cK.xyz (K = 0 is the lowest in xTB energy) and a CSV
for screen.py (names NAME_cK; the SMILES is only used if a geometry file is missing).
Usage: python conformers.py NAME [--n 80] [--window 3.0] [--max 4] [--out runs/conf]
"""
import argparse, csv, os, sys
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(here)); sys.path.insert(0, here)
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign
from psscreen.geometry import xtb_optimize, write_xyz
from molecules import SMILES

ap = argparse.ArgumentParser()
ap.add_argument("name"); ap.add_argument("--n", type=int, default=80)
ap.add_argument("--window", type=float, default=3.0, help="kcal/mol above the lowest xTB energy")
ap.add_argument("--max", type=int, default=4); ap.add_argument("--out", default=os.path.join(here, "runs", "conf"))
ap.add_argument("--rmsd", type=float, default=0.5, help="heavy-atom RMSD (A) below which conformers are duplicates")
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)

mol = Chem.AddHs(Chem.MolFromSmiles(SMILES[a.name]))
ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=a.n, randomSeed=11, pruneRmsThresh=0.3))
print(f"{a.name}: {len(ids)} embedded conformers", flush=True)
sym = [x.GetSymbol() for x in mol.GetAtoms()]
res = []
for i, cid in enumerate(ids):
    AllChem.MMFFOptimizeMolecule(mol, confId=cid, maxIters=1000)
    s, xyz, e = xtb_optimize(sym, mol.GetConformer(cid).GetPositions())
    res.append((e, xyz)); print(f"  conf {i}: E = {e:.4f} eV", flush=True)
res.sort(key=lambda t: t[0])
e0 = res[0][0]; EV_KCAL = 23.0605
heavy = [i for i, x in enumerate(sym) if x != "H"]
kept = []
for e, xyz in res:
    if (e - e0) * EV_KCAL > a.window or len(kept) >= a.max:
        continue
    # RMSD without alignment (rotation) is meaningless, so use Kabsch on heavy atoms
    def rmsd(p, q):
        p = p[heavy] - p[heavy].mean(0); q = q[heavy] - q[heavy].mean(0)
        u, _, vt = np.linalg.svd(p.T @ q); d = np.sign(np.linalg.det(u @ vt))
        r = u @ np.diag([1, 1, d]) @ vt
        return float(np.sqrt(((p @ r - q) ** 2).sum(1).mean()))
    if all(rmsd(xyz, k[1]) > a.rmsd for k in kept):
        kept.append((e, xyz))
rows = []
for k, (e, xyz) in enumerate(kept):
    nm = f"{a.name}_c{k}"
    write_xyz(os.path.join(a.out, f"{nm}.xyz"), sym, xyz, f"{nm} GFN2-xTB dE={(e - e0) * EV_KCAL:.2f} kcal/mol")
    rows.append((nm, SMILES[a.name]))
    print(f"kept {nm}: dE = {(e - e0) * EV_KCAL:.2f} kcal/mol", flush=True)
with open(os.path.join(a.out, f"{a.name}_conformers.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["name", "smiles"]); w.writerows(rows)
