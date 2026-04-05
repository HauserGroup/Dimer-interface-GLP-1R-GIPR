"""
"""
import pandas as pd 
import pathlib as pl 


def get_chain_from_residue_info(residue_str):
    """
    ['LEU $stru.A:279@O 307 \xa0\xa0\xa0',
    'TYR $stru.A:289@pi0 477 \xa0\xa0\xa0',
    'LEU $stru.A:279 305@CA \xa0\xa0\xa0',
    """
    # get chain (A)
    chain = residue_str.partition(".")[2].partition(":")[0]
    return chain


def get_residue_from_residue_info(residue_str):
    """
    return LEU 279
    """
    residue_3letter = residue_str.split()[0]
    residue_num = residue_str.split(":")[1].split()[0].split("@")[0]
    return f"{residue_3letter} {residue_num}"

if __name__ == "__main__":
    script_dir = pl.Path(__file__).parent.absolute()
    # Index(['Atom or residue 1', 'Atom or residue 2', 'Molecule Type', 'Interaction type', 'Distance', 'Distance type', 'Min Distance(Å)', 'Num Contacts']
    df = pd.read_excel(script_dir / "ic3nd_contacts.xlsx")
    df['Chain 1'] = df['Atom or residue 1'].apply(get_chain_from_residue_info)
    df['Chain 2'] = df['Atom or residue 2'].apply(get_chain_from_residue_info)
    df['Residue 1'] = df['Atom or residue 1'].apply(get_residue_from_residue_info)
    df['Residue 2'] = df['Atom or residue 2'].apply(get_residue_from_residue_info)
    
    # save Chain 1 residue numbers
    with open(script_dir / "ic3nd_contact_chainA_residues.txt", "w") as f:
        residues = df[df['Chain 1'] == 'A']['Residue 1'].unique()
        for r in residues:
            n = r.split()[1]
            f.write(f"{n}+")
    
    # save Chain 2 residue numbers
    with open(script_dir / "ic3nd_contact_chainB_residues.txt", "w") as f:
        residues = df[df['Chain 2'] == 'B']['Residue 2'].unique()
        for r in residues:
            n = r.split()[1]
            f.write(f"{n}+")