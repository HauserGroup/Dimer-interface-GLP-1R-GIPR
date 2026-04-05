"""
"""
import pathlib
from pymol import cmd 

MODELS = pathlib.Path("/Users/kcd635/HauserGroup Dropbox/luuk hd/GLP1R_GIPR_dimers/AF2-dimer-deepsampling/GLP1R_GIPR_dimers_originalseq")

for SEED in MODELS.iterdir():
    if not SEED.is_dir():
        continue
    for MODEL in SEED.iterdir():
        MODEL_NAME = MODEL.stem
        if MODEL_NAME.endswith('GLP-1'):
            PDB_P = MODEL / f"{MODEL_NAME}.pdb"
            print(f"Loading {PDB_P}")
            cmd.load(str(PDB_P), SEED)

