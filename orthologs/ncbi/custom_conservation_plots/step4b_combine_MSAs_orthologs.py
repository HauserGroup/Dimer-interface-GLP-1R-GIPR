"""
Combine GLP1R and GIPR MSAs into one file if species match.
"""
from step4_species_and_residues import read_fasta, get_id_to_seq
import pandas as pd
import pathlib as pl 


def sanitize_string(s : str):
    return s.replace(' ', '_').replace('/', '_').replace('-', '_').replace('.', '_').replace(',', '_').replace('(', '_').replace(')', '_').replace("'", '_').replace('"', '_').replace(';', '_').replace(':', '_').replace('?', '_').replace('!', '_')
                                                                                                                                                                                                                                                                                                                                             

if __name__ == "__main__":
    script_dir = pl.Path(__file__).parent.absolute()
    species = script_dir / 'step4a_species_and_ids.csv'
    species_df = pd.read_csv(species, index_col=0)
    # species,id_glp1r,id_gipr
    # Castor canadensis,XP_020041914.1,XP_020028513.1
    # Prionailurus viverrinus,XP_047716141.1,XP_047691825.1
    glp1r_aligned_msa = script_dir / 'GLP1R_clustalO/clustalo-I20250127-132055-0631-4283601-p1m/clustalo-I20250127-132055-0631-4283601-p1m.sqr.fa'
    glp1r_aligned_msa_records = read_fasta(glp1r_aligned_msa)
    glp1r_id_to_seq = get_id_to_seq(glp1r_aligned_msa_records)
    gipr_aligned_msa = script_dir / 'GIPR_clustalO/clustalo-I20250127-132518-0827-92000744-p1m/clustalo-I20250127-132518-0827-92000744-p1m.sqr.fa'
    gipr_aligned_msa_records = read_fasta(gipr_aligned_msa)
    gipr_id_to_seq = get_id_to_seq(gipr_aligned_msa_records)
    
    combined_msa = script_dir / 'step4b_combined_msa.fa'
    with open(combined_msa, 'w') as f:
        for i, rows in species_df.iterrows():
            glp1r_id = rows['id_glp1r']
            gipr_id = rows['id_gipr']
            species = i
            combined_id = f'{species}__{glp1r_id}__{gipr_id}'
            combined_id = sanitize_string(combined_id)
            f.write(f'>{combined_id}\n')
            glp1r_seq = glp1r_id_to_seq[glp1r_id]
            gipr_seq = gipr_id_to_seq[gipr_id]
            combined_seq = f'{glp1r_seq}{gipr_seq}'
            f.write(f'{combined_seq}\n')


    print(f'Combined MSA written to {combined_msa}')