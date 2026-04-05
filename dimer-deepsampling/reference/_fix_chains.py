from Bio.PDB import PDBParser, PDBIO
import pathlib


def update_chains(pdb_p, chain_remap, out_p):
    pdb_p = str(pdb_p)
    out_p = str(out_p)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_p)

    # Biopython structures are: Structure → Model → Chain → Residue → Atom
    # Create a copy of the structure with updated chains
    for model in structure:
        for chain in model:
            if chain.id in chain_remap:
                new_chain_id = chain_remap[chain.id]
                chain.id = new_chain_id
                

    io = PDBIO()
    io.set_structure(structure)
    io.save(out_p)


def main():
    """
    Original screen placed GLP1R -> GIPR
    Deepsampling places GIPR -> GLP1R

    Fix the reference, to be GIPR (B) -> GLP1R (C)
    """
    script_dir = pathlib.Path(__file__).parent.resolve()
    _fixed_chains = {
        'B' : 'B',
        'A' : 'C',
        'C' : 'D',
    }

    pdb_p = script_dir / 'full_w_glp1.pdb'    
    out_p = script_dir / 'full_w_glp1_updated_chains.pdb'
    update_chains(pdb_p, _fixed_chains, out_p)



if __name__ == "__main__":
    main()

