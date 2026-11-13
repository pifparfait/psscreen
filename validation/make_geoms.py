import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from psscreen.geometry import smiles_to_xyz, xtb_optimize, write_xyz
from molecules import SMILES
os.makedirs("runs", exist_ok=True)
for name, smi in SMILES.items():
    sym, xyz = smiles_to_xyz(smi)
    sym, xyz, e = xtb_optimize(sym, xyz)
    write_xyz(f"runs/{name}.xyz", sym, xyz, f"{name} GFN2-xTB {smi}")
    print(name, len(sym), e, flush=True)
