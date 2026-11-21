"""Chart 1: structures of the fifteen dyes, drawn with RDKit from the SMILES stored in the run files.

Usage (repository root): python paper/analysis/chart_structures.py
Output: paper/figures/chart_structures.png and .pdf
"""
import json
from pathlib import Path

from PIL import Image
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "validation" / "runs"
CAND = ROOT / "screening" / "runs" / "candidates_round1"
DYES = [("BODIPY_core", "BODIPY core"), ("BDP", "BDP"), ("I-BDP", "I-BDP"), ("I2-BDP", "I2-BDP"),
        ("I3-BDP", "I3-BDP"), ("azaBDP_Ph4", "aza-BDP-Ph4"), ("BDP-0", "BDP-0"), ("I-BDP-0", "I-BDP-0"),
        ("I2-BDP-0", "I2-BDP-0"), ("Cl2-BDP", "Cl2-BDP*"), ("Br2-BDP", "Br2-BDP*"), ("I1-BDP", "I1-BDP*"),
        ("azaBDP_Br2", "aza-BDP-Br2*"), ("azaBDP_BrPh2", "aza-BDP-BrPh2*"), ("azaBDP_IPh2", "aza-BDP-IPh2*")]


def smiles(name):
    for p in (RUNS / f"{name}.json", RUNS / "set3" / f"{name}_c0.json", CAND / f"{name}.json"):
        if p.exists():
            return json.loads(p.read_text())["smiles"]
    raise FileNotFoundError(name)


def main():
    mols, legends = [], []
    for name, label in DYES:
        m = Chem.MolFromSmiles(smiles(name))
        AllChem.Compute2DCoords(m)
        mols.append(m)
        legends.append(label)
    opts = Draw.MolDrawOptions()
    opts.legendFontSize = 30
    opts.bondLineWidth = 2
    opts.minFontSize = 16
    opts.addStereoAnnotation = False
    img = Draw.MolsToGridImage(mols, molsPerRow=5, subImgSize=(420, 340), legends=legends, drawOptions=opts)
    out = ROOT / "paper" / "figures" / "chart_structures"
    img.save(out.with_suffix(".png"))
    Image.open(out.with_suffix(".png")).convert("RGB").save(out.with_suffix(".pdf"), resolution=300)
    print("wrote", out.with_suffix(".png"))


if __name__ == "__main__":
    main()
