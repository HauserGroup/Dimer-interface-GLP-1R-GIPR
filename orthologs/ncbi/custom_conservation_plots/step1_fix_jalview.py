"""
1. From the jalview export, only take the residues of GLP1R.
"""
import pandas as pd
import pathlib as pl
import seaborn as sns
import matplotlib.pyplot as plt
sns.set_style("whitegrid")


def read_jalview_annotations_export(p):
    """
    in jalview, go to file -> export -> annotations and save as .csv
    """
    annot = pd.read_csv(p, index_col=0, delimiter=',', header=None)
    annot = annot.T
    # if two times 'Consensus' in the header, rename second
    duplicated_cols = annot.columns.duplicated()
    annot.columns = [f'{col} {i}' if duplicated else col for i, (col, duplicated) in enumerate(zip(annot.columns, duplicated_cols))]
    # strip each col and val
    annot.columns = [col.strip() for col in annot.columns]
    # detect nan as empty string (after strip)
    annot = annot.applymap(str.strip)
    # set '' to nan
    annot = annot.replace('', pd.NA)
    nan_rows = annot.isna().all(axis=1)
    annot = annot[~nan_rows]
    # make numeric if possible
    annot = annot.apply(pd.to_numeric, errors='ignore')
    # save as fixed
    #annot.to_csv(p.with_suffix('.fixed.csv'))
    return annot


def add_amino_acids_per_row(protein_seq, amino_acid_col, df):
    """
    """
    assert amino_acid_col not in df.columns, 'amino_acid_col already in df'
    prot_len = len(protein_seq)
    rows = df.shape[0]
    assert prot_len == rows, f"protein_seq length {prot_len} != df rows {rows}"

    df[amino_acid_col] = list(protein_seq)
    return df


#>NP_002053.3
glp1r="""
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
gipr="""
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
    # settings
    glp1r_chains = pd.read_csv('/Users/kcd635/HauserGroup Dropbox/luuk hd/GLP1R_GIPR_dimers/code/glp1r_human_chain_info.csv')
    glp1r_tm_residues = glp1r_chains[glp1r_chains['protein_segment'].str.startswith('TM')]['sequence_number'].values.tolist()

    gipr_chains = pd.read_csv('/Users/kcd635/HauserGroup Dropbox/luuk hd/GLP1R_GIPR_dimers/code/gipr_human_chain_info.csv')
    gipr_tm_residues = gipr_chains[gipr_chains['protein_segment'].str.startswith('TM')]['sequence_number'].values.tolist()

    # 279+289+282+283+286+287+290+291+275 
    glp1r_construct_residues = [279,289,282,283,286,287,290,291,275]
    glp1r_interface_residues = [279,289,282,283,286,287,290,291]
    # 136+147+192+150+151+146+142+143+188+139+191+138+141 
    gipr_construct_residues = [136,147,192,150,151,146,142,143,188,139,191,138,141]
    gipr_interface_residues = [136,147,192,150,151,146,142,143,188,139,191,138,141]

    # fix glp1r
    script_dir = pl.Path(__file__).parent
    annot = read_jalview_annotations_export(script_dir / "GLP1R_jalview_annotation.csv")
    annot = add_amino_acids_per_row(glp1r, 'Residue', annot)
    annot = annot[annot['Residue'] != '-']
    annot = annot.reset_index(drop=True)
    annot['sequence_number'] = annot.index + 1
    annot['in_TM'] = annot['sequence_number'].isin(glp1r_tm_residues)
    annot['in_construct'] = annot['sequence_number'].isin(glp1r_construct_residues)
    annot['in_interface'] = annot['sequence_number'].isin(glp1r_interface_residues)
    annot.to_csv(script_dir / "GLP1R_jalview_annotation_fixed.csv")

    # fix GIPR
    annot = read_jalview_annotations_export(script_dir / "GIPR_jalview_annotation.csv")
    annot = add_amino_acids_per_row(gipr, 'Residue', annot)
    annot = annot[annot['Residue'] != '-']
    annot = annot.reset_index(drop=True)
    annot['sequence_number'] = annot.index + 1
    annot['in_TM'] = annot['sequence_number'].isin(gipr_tm_residues)
    annot['in_construct'] = annot['sequence_number'].isin(gipr_construct_residues)
    annot['in_interface'] = annot['sequence_number'].isin(gipr_interface_residues)
    annot.to_csv(script_dir / "GIPR_jalview_annotation_fixed.csv")

