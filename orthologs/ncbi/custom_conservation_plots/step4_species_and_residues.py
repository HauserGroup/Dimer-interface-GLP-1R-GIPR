"""
1. match proteins by species.

"""
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pathlib
from Bio import SeqIO


def read_fasta(fasta_file):
    """ Read fasta file and return a dictionary with the sequence.
    """
    records = []
    for record in SeqIO.parse(fasta_file, 'fasta'):
        records.append(record)
    return records


def get_species_from_description(description):
    """ Extract species from the description.
    """
    return description.split('[')[-1].split(']')[0]


def get_id_to_species(ortho_records):
    """ Extract species from the id.
    
    [
    {'_seq': Seq('MTTSPILQLLLRLSLCGLLLQRAETGSKGQTAGELYQRWERYRRECQETLAAAE...SYC'), 
     'id': 'NP_000155.1', 
     'name': 'NP_000155.1', 
     'description': 'NP_000155.1 gastric inhibitory polypeptide receptor isoform 1 precursor [Homo sapiens]',
     'dbxrefs': [],
     'annotations': {},
     '_per_letter_annotations': {},
     'features': []}
    ]
    """
    id_to_species = dict()
    for record in ortho_records:
        description = record.description
        species = get_species_from_description(description)
        assert record.id not in id_to_species, 'Already parsed'
        id_to_species[record.id] = species
    
    return id_to_species


def get_id_to_seq(records):
    id_to_seq = dict()
    for r in records:
        id_to_seq[r.id] = str(r.seq)
    return id_to_seq


def fill_df_based_on_alignments(df,
                                aligned_query,
                                id_to_species,
                                alignments,
                                protein_id,
    ):
    ""
    for id, species in id_to_species.items():
        aligned_target = alignments[id]
        query_target_zipped = zip(aligned_query, aligned_target)
        i = 0
        for q, t in query_target_zipped:
            if q == '-':
                continue
            i += 1
            col = f'{protein_id}_{q}{i}'
            assert col in df.columns, f'Column {col} not in df'
            df.loc[df['species'] == species, col] = t
    return df


#>NP_002053.3
glp1r_aln_seq="""
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
-------------------------MA----GAP------GPLRLAL-------------
------------LL---LGMVGR--------AGPR-PQ----------------------
----------GATVSLWETVQKWREYRRQCQRS--------LTEDPPPAT----------
-------DLFCNRTFDEYACWPDGEPGSF----V----------NVSCPWYLPWAS----
-SVPQG--HVYRFCTAEGLWLQK---DNSSLPWRDLSECEESKR-GERSSPEEQLLFLYI
IYTVGYALSFSALVIASAILLGFR------------------HLHCT---RNYIHLNLFA
SFILRALSVFIKDAALKWMYSTAAQQHQWDGLLSYQDSLSCRLVFLLMQYCVAANYYWLL
VEGVYLYTLLAFSVLSEQWIFRLYVSIGWGVPLLFVVPWGIVKYLYEDEGCWTRNSNMNY
WLIIRLPILFAIGVNFLIFVRVICIVVSKLKANLMCKTDIKCR-----------------
--LAKSTLTLIPLLGTHEVIFAFVMDEHARGTLRFIKLFTELSFTSF-------------
------------------------------------------------------------
------------------QGLMVAILYCFVNNE-VQLEF---------RKSWER-WRLEH
LHI---------QRD---SS-------MKPLK----CP--------------------TS
SLSSGAT-AGSSM----------YTAT-CQASCS--------------------------
------------------------------------------------------------
-------------------
""".replace("\n","").strip()


# >NP_000155.1
gipr_aln_seq="""
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
------------------------------------------------------------
--------------------MTTSPILQLL----LRLSLCGLLLQRA-------------
ETGSKGQTAGE-LYQRWERYRRECQETLAAAEPP--SG----LACNGSFDMYVCWDYAAP
NATARASCPWYLPWHHHV------------------------------------AAGFVL
RQCGSDGQWG------L-----------------------------W-----RDHTQCEN
PEKNEAFL------------------------------------DQRLILERLQVMYTVG
YSLSLATLLLALLILSLFRRLHCTR-NYIHINLFTSFMLRA--------------AAIL-
SRDRLLPRPGP----------YLGDQALALW--------------NQALAACRTAQI---
---------------VTQYCVGANYTWLLVEGVYLH--SLLVLVGGSE--------EGHF
RYYLLLGWGA--------------------------------------------------
---------------------P-ALFVIPWVIVRYLYE---NTQCWERNEVK--------
---------------------------------------AIWWII---------------
RTPILMTIL------------------------------------------INFLIFIRI
LGILLSKLRTRQMRCRDYRLRL--------ARSTLTLVPLLGVHEVVFAPV-TEEQA-RG
ALRFAKLGFEIFLSSFQ-------------GFLVSVL------YCFINK-----------
----------E-------------------------------------------------
---------------VQSEIRRG---WHHCRLRRSLGEEQR-QL---------P------
---------ERAF--------RAL------------------------------PSGSGP
G------------------------E------------------VPTSRGLSS--GTLPG
PGNEA-SRELESYC----------------------------------------------
------------------------------------------------------------
-------------------------------------------------
""".replace("\n","").strip()


if __name__ == "__main__":
    """
    """
    # get the species
    script_dir = pathlib.Path(__file__).parent
    ortho_glp1r = script_dir / 'GLP1R_refseq_protein.fasta'
    ortho_gipr = script_dir / 'GIPR_refseq_protein.fasta'
    glp1r_id_to_species = get_id_to_species(read_fasta(ortho_glp1r))
    glp1r_species = set(glp1r_id_to_species.values())
    gipr_id_to_species = get_id_to_species(read_fasta(ortho_gipr))
    gipr_species = set(gipr_id_to_species.values())

    # get the overlap we can keep
    overlap = glp1r_species.intersection(gipr_species)
    # filter gipr_id_to_species and glp1r_id_to_species
    glp1r_id_to_species = {k: v for k, v in glp1r_id_to_species.items() if v in overlap}
    gipr_id_to_species = {k: v for k, v in gipr_id_to_species.items() if v in overlap}

    # create a dataframe with [species, glp1r id, gipr id]
    df = pd.DataFrame(columns=['species', 'id_glp1r', 'id_gipr'])
    for species in overlap:
        glp1r_id = [k for k, v in glp1r_id_to_species.items() if v == species]
        gipr_id = [k for k, v in gipr_id_to_species.items() if v == species]
        species_df = pd.DataFrame({'species': species, 'id_glp1r': glp1r_id, 'id_gipr': gipr_id})
        df = pd.concat([df, species_df], axis=0)
    df.to_csv(script_dir / 'step4a_species_and_ids.csv', index=False)

    # add the residues from the alignments
    gipr_alignments = script_dir / 'GIPR_clustalO/clustalo-I20250127-132518-0827-92000744-p1m/clustalo-I20250127-132518-0827-92000744-p1m.sqr.fa'
    gipr_alignments = get_id_to_seq(read_fasta(gipr_alignments))
    glp1r_alignments = script_dir / 'GLP1R_clustalO/clustalo-I20250127-132055-0631-4283601-p1m/clustalo-I20250127-132055-0631-4283601-p1m.sqr.fa'
    glp1r_alignments = get_id_to_seq(read_fasta(glp1r_alignments))
    # use aln seq to get the residues that are not - in the alignment
    gipr_indexes_to_keep = [i for i, aa in enumerate(gipr_aln_seq) if aa != '-']
    glp1r_indexes_to_keep = [i for i, aa in enumerate(glp1r_aln_seq) if aa != '-']

    # create a df with the residues there. Species as index. 
    # with columns as 'glp1r_M1', 'glp1r_A2', ..., 'gipr_M1', 'gipr_T2'
    gipr_seq = gipr_aln_seq.replace('-', '')
    glp1r_seq = glp1r_aln_seq.replace('-', '')
    # add cols
    new_cols =  [f'gipr_{aa}{i+1}' for i, aa in enumerate(gipr_seq)]
    new_cols = new_cols + [f'glp1r_{aa}{i+1}' for i, aa in enumerate(glp1r_seq)]
    df = pd.concat([df, pd.DataFrame(columns=new_cols)], axis=1)
    df = fill_df_based_on_alignments(df.copy(), gipr_aln_seq, gipr_id_to_species, gipr_alignments, 'gipr')
    df = fill_df_based_on_alignments(df.copy(), glp1r_aln_seq, glp1r_id_to_species, glp1r_alignments, 'glp1r')
    df.to_csv(script_dir / 'step4b_species_and_residues.csv', index=False)

