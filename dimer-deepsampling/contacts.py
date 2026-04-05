"""
From a .cif or .pdb file, get the contacts (residue distances in angstrom).
"""
import pathlib
import pandas as pd
import numpy as np
import itertools
import pathlib as pl
import pandas as pd
import numpy as np 
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import protein_letters_3to1 as three_to_one
from Bio.PDB import PDBParser
import re
from io import StringIO


def load_structure(pdb_path: pl.Path):
    """
    Load a PDB structure after preprocessing its content in memory.
    Cleans any non-integer residue IDs or insertion codes on the fly.
    """
    cleaned_lines = []
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Strip non-integer characters (e.g., insertion codes) from residue sequence ID
                line = re.sub(r"(\s\d+)[A-Z]?", r"\1", line)
            cleaned_lines.append(line)
    
    # Join the cleaned lines and parse as a file-like object
    cleaned_pdb_content = StringIO("".join(cleaned_lines))
    
    # Use Biopython's PDBParser to load the cleaned content
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', cleaned_pdb_content)
    return structure


def get_chain_to_sequence(model):
    """
    Get the sequence of the chain in the model.
    Handles residues with non-integer sequence IDs and insertion codes.
    
    Parameters:
    model: Biopython Model object
    
    Returns:
    dict: A dictionary mapping chain IDs to their respective sequences.
    """
    chain_to_seq = {}
    for chain in model:
        sequence = []
        for residue in chain.get_residues():
            # Skip non-standard residues (e.g., water or HETATM)
            if residue.id[0] != ' ':
                continue
            try:
                resname = residue.get_resname()
                # Convert three-letter code to one-letter code
                one_letter = three_to_one[resname]
                sequence.append(one_letter)
            except KeyError:
                # Handle unknown residues or non-standard amino acids
                sequence.append('X')
            except ValueError as e:
                # print(f"Warning: Skipping residue {residue} due to error: {e}")
                continue
        # Add sequence for this chain
        chain_to_seq[chain.id] = ''.join(sequence)
    return chain_to_seq


def get_residues_for_chain(loaded_struct, chain_letter : str):
    """
    Get residues for a chain letter.
    Returns a list of residues.
    """
    for chain in loaded_struct.get_chains():
        if chain.id == chain_letter:
            return [residue for residue in chain.get_residues()]
    
    raise ValueError(f"Chain {chain_letter} not found in structure {loaded_struct.id}")


def get_residue_names_for_chain(loaded_struct, chain_letter : str):
    """
    Get residue numbers for a chain letter.
    Returns a list of residue numbers.
    """
    residues = get_residues_for_chain(loaded_struct, chain_letter)
    residue_numbers = [residue.id[1] for residue in residues]
    residue_names = [residue.get_resname() for residue in residues]
    residue_names = list(map(lambda x: three_to_one[x], residue_names))
    return [f"{name}{num}" for name, num in zip(residue_names, residue_numbers)]


def get_alpha_carbon_cmap(loaded_struct, chain_a, chain_b):
    """
    return cmap of angstrom distances 
    """
    # Get residues 
    chain_a_residues = get_residues_for_chain(loaded_struct, chain_a)
    chain_b_residues = get_residues_for_chain(loaded_struct, chain_b)
    # Get residue numbers
    chain_a_resnames = get_residue_names_for_chain(loaded_struct, chain_a)
    chain_b_resnames = get_residue_names_for_chain(loaded_struct, chain_b)
    # Get alpha-carbon atoms
    chain_a_ca_atoms = [residue['CA'] for residue in chain_a_residues if 'CA' in residue]
    chain_b_ca_atoms = [residue['CA'] for residue in chain_b_residues if 'CA' in residue]
    # distance cmap using numpy
    chain_a_ca_coords = np.array([atom.coord for atom in chain_a_ca_atoms])
    chain_b_ca_coords = np.array([atom.coord for atom in chain_b_ca_atoms])
    if chain_b_ca_coords.size == 0 or chain_a_ca_coords.size == 0:
        raise ValueError(f"No CA atoms found in chain {chain_b} or chain {chain_a}")
    # euclidean distance between all pairs of CA atoms
    cmap = np.linalg.norm(chain_b_ca_coords[:, np.newaxis] - chain_a_ca_coords, axis=-1)
    if cmap.size == 0:
        raise ValueError(f"Distance cmap is empty for chain {chain_b} and chain {chain_a}")
    # transform to df, with residue names
    cmap = pd.DataFrame(cmap, index=chain_b_resnames, columns=chain_a_resnames)
    return cmap


def _get_chain_distance_maps(structure, distance_method="CA"):
    """
    return dict, with tuple to pd DataFrame distance maps between chains
    """
    chain_to_seq = get_chain_to_sequence(structure)
    chains = chain_to_seq.keys()
    ###print(f"Chains found in structure: {chains}")
    if len(chains) < 2:
        raise ValueError(f"Structure must contain at least two chains to compute inter-chain distances. Found {len(chains)} chain(s).")

    # get chain combinations (unique, sorted)
    chain_combinations = list(itertools.combinations(sorted(chains), 2))
    ##print(f"Chain combinations to analyze: {chain_combinations}")

    # get cmap for each chain combination
    out = {}
    for chain_a, chain_b in chain_combinations:
        if distance_method == "CA":
            cmap = get_alpha_carbon_cmap(structure, chain_a, chain_b)
            # report any nan
            if cmap.isnull().values.any():
                raise ValueError(f"NaN values found in distance map between chain {chain_a} and chain {chain_b}")
            out[(chain_a, chain_b)] = cmap
        else:
            raise ValueError(f"Unknown distance method: {distance_method}")

    return out


def get_angstrom_distances_for_pdb(pdb_path: pl.Path, distance_method="CA"):
    """
    Get angstrom distances for a pdb file.
    Returns pd DataFrame distance maps for each chain pair in the structure.
    """
    structure = load_structure(pdb_path)
    structure = structure[0]  # get first model
    distances = _get_chain_distance_maps(structure, distance_method=distance_method)
    return distances


def example():
    script_dir = pathlib.Path(__file__).parent
    # AF2-dimer-deepsampling/GLP1R_GIPR_dimers/
    pdb_p = script_dir / 'GLP1R_GIPR_dimers/seed1/gipr_human___glp1r_human___Exendin-4/ranked_0.pdb'   
    if not pdb_p.exists():
        raise FileNotFoundError(f"Example PDB file not found at {pdb_p}. Please ensure the file exists.")
    
    structure = load_structure(pdb_p)
    structure = structure[0]
    distances = _get_chain_distance_maps(structure, distance_method="CA")
    print(distances)


def ref():
    script_dir = pathlib.Path(__file__).parent
    pdb_p = script_dir / 'reference/full_w_glp1.pdb'
    structure = load_structure(pdb_p)
    structure = structure[0]
    distances = _get_chain_distance_maps(structure, distance_method="CA")
    # save to csv
    for (chain_a, chain_b), cmap in distances.items():
        out_p = script_dir / f'reference/reference_distance_map_{chain_a}_{chain_b}.csv'
        cmap.to_csv(out_p)
        print(f"Saved distance map for chains {chain_a} and {chain_b} to {out_p}")


if __name__ == "__main__":
    print("Running b_contacts.py example...")
    # example()
    ref()
