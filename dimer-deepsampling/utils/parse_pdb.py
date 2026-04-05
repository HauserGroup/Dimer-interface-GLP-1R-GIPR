""" For each structure, 
Get the peptide length (resolved)
Family.
Class.
"""
import pathlib as pl
import pandas as pd
import seaborn as sns 
import matplotlib.pyplot as plt
import numpy as np 
from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.Polypeptide import protein_letters_3to1 as three_to_one

# add parent
import sys
sys.path.append(str(pl.Path(__file__).parent.parent.absolute()))

from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB import PDBParser
import re
from io import StringIO

# pairwise alignment
from Bio import pairwise2
from Bio.pairwise2 import format_alignment



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



def get_lengths_of_pdbs(pdb_dir: pl.Path):
    script_dir = pl.Path(__file__).parent
    structs = pdb_dir.glob("*.pdb")
    out = []
    for s in structs:
        name = s.stem
        # print(name)
        model = load_structure(s)
    
        chain_to_seq = get_chain_to_sequence(model)
        assert len(chain_to_seq) == 1, f"!= 1 chain {s} {chain_to_seq}"
        seq = list(chain_to_seq.values())[0]
        l = len(seq)
        out.append((name, l, seq))
        
    df = pd.DataFrame(out, columns=['peptide', 'l', 'sequence'])
    return df

