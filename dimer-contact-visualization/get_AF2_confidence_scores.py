"""
In the ./raw/ dir there are subdirs with the pdb files and ranking_debug.json files
from AF predictions.

This code will look at the ranking scores and create a plot of how GPCR dimers differ in iptm+ptm
scores. 
"""
import pathlib
import os
import json
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from Bio import PDB


def get_iptm_ptm_scores(path_to_ranking_json):
    """
    example returns:
    {'model_1_multimer_v3_pred_0': 0.3314535102120189, 'model_2_multimer_v3_pred_0': 0.3177772733702266, 'model_3_multimer_v3_pred_0': 0.43851442414153263, 'model_4_multimer_v3_pred_0': 0.2625854925756743, 'model_5_multimer_v3_pred_0': 0.3936904594183363}
    """
    with open(path_to_ranking_json, "r") as f:
        ranking = json.load(f)
    return ranking["iptm+ptm"]


def get_highest_iptm_ptm_score(iptm_ptm_scores):
    """
    example returns:
    ('model_3_multimer_v3_pred_0', 0.43851442414153263)
    """
    id_of_best, score_of_best = max(iptm_ptm_scores.items(), key=lambda x: x[1])
    return id_of_best, score_of_best


def get_iptm_scores_from_AF_dirs(af_dirs):
    # create empty df with columns: gpcr_a, gpcr_b, peptide, id_of_best, iptm_ptm_score_of_best
    df = pd.DataFrame(columns=["gpcr_a", "gpcr_b", "peptide", "id_of_best", "iptm_ptm_score_of_best"])
    
    # get the ranking scores
    for af_dir in af_dirs:
        ranking_json = list(af_dir.glob("ranking_debug.json"))[0]
        ranking_scores = get_iptm_ptm_scores(ranking_json)
        highest_id, highest_score = get_highest_iptm_ptm_score(ranking_scores)

        parts = af_dir.name.split("___")
        peptide = None
        if len(parts) == 2:
            gpcr_a = parts[0]
            gpcr_b = parts[1]
        elif len(parts) == 3:
            gpcr_a = parts[0] 
            gpcr_b = parts[1]
            peptide = parts[2]
            if peptide not in ["gip", "glp1"]:
                peptide = peptide.split("_")[0]
        else:
            raise ValueError(f"AF_DIR.name {af_dir.name} has wrong number of parts")
        
        # save to df
        df = df._append({"gpcr_a": gpcr_a, "gpcr_b": gpcr_b, "peptide": peptide, "id_of_best": highest_id, "iptm_ptm_score_of_best": highest_score}, ignore_index=True)
    return df


def calculate_overall_pLDDT(pdb_path):
    """
    """
    # get the pLDDT scores
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("pdb", pdb_path)
    model = structure[0]
    pLDDTs = []
    for chain in model:
        for residue in chain:
            for atom in residue:
                if atom.get_name() == "CA":
                    pLDDTs.append(atom.bfactor)
    # calculate the mean
    return np.mean(pLDDTs)


def get_plddt_df(pdb_list):
    """
    """
    df = pd.DataFrame(columns=["gpcr_a", "gpcr_b", "peptide", "pLDDT"])
    for pdb in pdb_list:
        name_stripped = pdb.name.split(".")[0]
        parts = name_stripped.split("___")
        gpcr_a = parts[0]
        gpcr_b = parts[1]
        if len(parts) > 2:
            peptide = parts[2].split("_")[0]
        else:
            peptide = np.NaN
        pLDDT = calculate_overall_pLDDT(pdb)
        df = df._append({"gpcr_a": gpcr_a,
                        "gpcr_b": gpcr_b,
                        "peptide": peptide, 
                        "pLDDT": pLDDT}, ignore_index=True)
    return df


if __name__ == "__main__":
    RAW_DIR = pathlib.Path("./data/AF2_predictions_raw/")
    PDB_DIR = pathlib.Path("./data/AF2_PDBs/")
    # Get all AF2 prediction dirs, in this case containing "gipr_human" or "glp1r_human"
    AF_DIRS = [x for x in RAW_DIR.iterdir() if x.is_dir() and "gipr_human" in x.name or "glp1r_human" in x.name]

    # Part 1, copy the PDBs to a new dir
    if not PDB_DIR.exists():
        PDB_DIR.mkdir()
    for AF_DIR in AF_DIRS:
        # the PDB file has the same name as the dir, but suffixed with .pdb
        PDB_PATH = AF_DIR / (AF_DIR.name + ".pdb")
        if os.path.exists(PDB_PATH):
            NEW_PDB_P = PDB_DIR / PDB_PATH.name
            # copy the file to the new dir
            os.system(f"cp {PDB_PATH} {NEW_PDB_P}")
        else:
            raise ValueError(f"AF_DIR {AF_DIR} has no PDB file", PDB_PATH)
        
    # Part 2, get the confidence score stats
    IPTM_STATS = pathlib.Path("./data/dimer_iptm_ptm_scores.csv")
    PLDDT_STATS = pathlib.Path("./data/dimer_pLDDT_scores.csv")
    if not IPTM_STATS.exists():
        df = get_iptm_scores_from_AF_dirs(AF_DIRS)
        df.to_csv(IPTM_STATS)

    # Part 3, get the pLDDT scores
    PDB_LIST = [x for x in PDB_DIR.iterdir() if x.is_file() and x.suffix == ".pdb" and "gipr_human" in x.name or "glp1r_human" in x.name]
    if not PLDDT_STATS.exists():
        df = get_plddt_df(PDB_LIST)
        df.to_csv(PLDDT_STATS)

