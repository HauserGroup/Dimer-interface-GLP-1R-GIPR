"""If AlphaFold2 has been trained on a receptor complex, then the complex is considered known.
Especially if the true complex with the (principal agonist) exists in the PDB.
This script looks for known complexes for the dataset.
"""

import requests
import json

import pandas as pd


def get_structures_data_from_GPCRdb(gpcr_id):
    """
    gpcr_id: str
        The GPCRdb ID for the GPCR. For example, 'oprm_human'
    """
    # get a list of structures of a protein
    # /services/structure/protein/{entry_name}/
    url = f"https://gpcrdb.org/services/structure/protein/{gpcr_id}/"
    response = requests.get(url)
    response.raise_for_status()
    data = response.json()
    return data


def get_structure_info_for_gpcr(gpcr_id):
    data = get_structures_data_from_GPCRdb(gpcr_id)
    df = pd.DataFrame(data)
    df.to_csv(f"{gpcr_id}_structures.csv", index=False)
    # publication_date -> 2022-11-09

    # ligands -> "[{'name': 'fentanyl', 'type': 'small-molecule', 'function': 'Agonist', 'PDB': '7V7', 'SMILES': 'CCC(=O)N(c1ccccc1)C1CCN(CC1)CCc1ccccc1'}]"


def test_get_GtP_ligands():
    # /services/ligands/gtp_ids/
    url = "https://gpcrdb.org/services/ligands/gtp_ids/"
    response = requests.get(url)
    response.raise_for_status()
    data = response.json()

    with open("test.tsv", "w") as f:
        for lig_dict in data:
            f.write(f"{lig_dict['gtp_ligand_id']} {lig_dict['ligand_name']}\n")


def get_chain_info(gpcr_id):
    """from GPCRdb, find out residues and chains of the GPCR

    returns: list of dicts
    where every Dict is info for 1 residue, like:
    {'sequence_number': 465, 'amino_acid': 'Y',
    'protein_segment': 'C-term', 'display_generic_number': None}
    """

    url = f"https://gpcrdb.org/services/residues/{gpcr_id}/"
    response = requests.get(url)
    response.raise_for_status()
    data = response.json()
    return data


def save_chain_info(gpcr_id, chain_path):
    data = get_chain_info(gpcr_id)
    df = pd.DataFrame(data)
    df.to_csv(chain_path, index=False)


if __name__ == "__main__":
    """
    """
    get_structure_info_for_gpcr("glp1r_human")
    get_structure_info_for_gpcr("gipr_human")

    save_chain_info("glp1r_human", "glp1r_human_chain_info.csv")
    save_chain_info("gipr_human", "gipr_human_chain_info.csv")
