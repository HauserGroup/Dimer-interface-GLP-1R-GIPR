"""
Parse all dimer predictions to a .csv
"""
import pandas as pd
import numpy as np
import pathlib
import logging
import json
import tqdm
import itertools
from concurrent.futures import ThreadPoolExecutor, as_completed

from contacts import get_angstrom_distances_for_pdb
MODEL_BASE_DIR = pathlib.Path(__file__).parent / 'GLP1R_GIPR_dimers_originalseq'
GLP1R_CHAIN = 'B'
# GLP1R_LENGTH = 418
GLP1R_LENGTH = 463 # full length including extra residues

GIPR_CHAIN = 'C'
# GIPR_LENGTH = 410
GIPR_LENGTH = 466 # full length including extra residues
CROP = False

def json_to_dict(json_p : pathlib.Path) -> dict:
    with open(json_p, 'r') as f:
        d = json.load(f)
    return d


def _parse_pdb_p(pdb_p : str) -> dict:
    """
    GLP1R_GIPR_dimers/seed3/gipr_human___glp1r_human___Exendin-9/unrelaxed_model_4_multimer_v3_pred_0.pdb
    
    seed: 3

    peptide: Exendin-9

    model: 4
    sample (prediction run): 0
    """
    pdb_p = pathlib.Path(pdb_p)

    # gipr_human___glp1r_human___Exendin-9
    parent_name = pdb_p.parent.name
    peptide = parent_name.split('___')[-1]

    # unrelaxed_model_4_multimer_v3_pred_0
    file_stem = pdb_p.stem
    sample = file_stem.split('_')[-1]
    model = file_stem.partition('_model_')[2].split('_')[0]

    # seed3
    seed = pdb_p.parent.parent.name.replace('seed', '')

    # get iptm
    expected_confidence_json = pdb_p.parent / 'ranking_debug.json'
    if not expected_confidence_json.exists():
        logging.error(f'Expected confidence json does not exist: {expected_confidence_json}')
        iptm = float('nan')
    else:
        confidence = json_to_dict(expected_confidence_json)
        confidence_key = file_stem.replace('unrelaxed_', '')
        iptm = confidence['iptm+ptm'][confidence_key]
    
    return {
        'stem': pdb_p.stem,
        'seed': int(seed),
        'peptide': peptide,
        'model': int(model),
        'sample': int(sample),
        'iptm': float(iptm),
    }
    

def _add_parsed_to_df(row):
    parsed = _parse_pdb_p(row['pdb_p'])
    for k, v in parsed.items():
        row[k] = v
    return row


def _get_id(row):
    # seed3/gipr_human___glp1r_human___Exendin-9/unrelaxed_model_4_multimer_v3_pred_0
    # -> seed3|Exendin-9|unrelaxed_model_4_multimer_v3_pred_0
    return f'seed{row["seed"]}|{row["peptide"]}|{row["stem"]}'    


def _process_angstrom_distances_for_row(row, residue_dist_dir: pathlib.Path):
    """
    Per chain combination, save a csv with distances.
    Use id + chainA + chainB to name file.
    """
    pdb_p = pathlib.Path(MODEL_BASE_DIR) / row['pdb_p']
    distances = get_angstrom_distances_for_pdb(pdb_p, distance_method="CA")

    # per chain combination
    for (chain_a, chain_b), group in distances.items():
        out_p = residue_dist_dir / f"{row['id']}_{chain_a}|{chain_b}.csv"
        logging.debug(f'Group shape for {chain_a}|{chain_b}: {group.shape}')
        group.to_csv(out_p, index=True)


def assign_receptor_bound_from_paths(gipr_csv: str, glp1r_csv: str) -> str:
    gipr = pd.read_csv(gipr_csv, index_col=0)
    gipr_dimension = gipr.shape.index(max(gipr.shape))
    gipr_residue_contacts = (gipr.min(axis=gipr_dimension) < 10.0).sum()

    glp1r = pd.read_csv(glp1r_csv, index_col=0)
    glp1r_dimension = glp1r.shape.index(max(glp1r.shape))
    glp1r_residue_contacts = (glp1r.min(axis=glp1r_dimension) < 10.0).sum()

    contact_diff = abs(int(gipr_residue_contacts) - int(glp1r_residue_contacts))
    if contact_diff < 10:
        return "Ambiguous"
    return "GIPR" if gipr_residue_contacts > glp1r_residue_contacts else "GLP1R"


def _assign_receptor_task(args):
    # top-level function so it’s picklable for multiprocessing
    idx, gipr_csv, glp1r_csv = args
    rb = assign_receptor_bound_from_paths(gipr_csv, glp1r_csv)
    return idx, rb


def prepare_distances(out_p: pathlib.Path, residue_dist_dir: pathlib.Path, pdb_d : pathlib.Path, _debug : bool = False):
    # load
    logging.info(f"Parsing PDB files in directory: {pdb_d}")
    pdbs = list(pdb_d.glob('*/*/*.pdb'))
    if _debug:
        # take first 10
        logging.warning("Debug mode: only processing first 10 PDB files.")
        pdbs = pdbs[:10]
    logging.info(f"Found {len(pdbs)} PDB files.")
    
    # create df for pdb_p, use relative paths
    df = pd.DataFrame({'pdb_p': [str(pdb) for pdb in pdbs]})
    # drop those without 'unrelaxed'
    df = df[df['pdb_p'].str.contains('unrelaxed')]
    df = df.reset_index(drop=True)
    df = df.apply(_add_parsed_to_df, axis=1)
    df['id'] = df.apply(_get_id, axis=1)
    assert df['id'].is_unique, "IDs are not unique!"


    # calculate angstrom distances in parallel
    logging.info("Calculating angstrom distances for each PDB...")
    with ThreadPoolExecutor(max_workers=8) as executor:
        futures = []
        for _, row in df.iterrows():
            futures.append(executor.submit(_process_angstrom_distances_for_row, row, residue_dist_dir))
        for future in tqdm.tqdm(as_completed(futures), total=len(futures)):
            try:
                future.result()
            except Exception as e:
                logging.error(f"Error processing a PDB: {e}")

    # move pdb_p to end
    df = df[[col for col in df.columns if col != 'pdb_p'] + ['pdb_p']]
    df['dimer_distances'] = str(residue_dist_dir) + '/seed' + df['seed'].astype(str) + '|' + df['peptide'] + '|' + df['stem'] + '_B|C.csv'
    df['GIPR_pep_distances'] = str(residue_dist_dir) + '/seed' + df['seed'].astype(str) + '|' + df['peptide'] + '|' + df['stem'] + f'_{GIPR_CHAIN}|D.csv'
    df['GLP1R_pep_distances'] = str(residue_dist_dir) +'/seed' + df['seed'].astype(str) + '|' + df['peptide'] + '|' + df['stem'] + f'_{GLP1R_CHAIN}|D.csv'
    df["receptor_bound"] = pd.NA

    # process with _assign_which_receptor_is_bound, use ThreadPoolExecutor
    logging.info("Assigning which receptor is bound for each PDB...")
    tasks = [
        (idx, row["GIPR_pep_distances"], row["GLP1R_pep_distances"])
        for idx, row in df[["GIPR_pep_distances", "GLP1R_pep_distances"]].iterrows()
    ]
    logging.info(f"Processing {len(tasks)} tasks for receptor bound assignment...")
    with ThreadPoolExecutor(max_workers=8) as executor:
        futures = {executor.submit(_assign_receptor_task, task): task[0] for task in tasks}
        for future in tqdm.tqdm(as_completed(futures), total=len(futures)):
            idx = futures[future]
            try:
                result_idx, rb = future.result()
                df.at[result_idx, 'receptor_bound'] = rb
            except Exception as e:
                logging.error(f"Error assigning receptor bound for index {idx}: {e}")

    # save
    df.to_csv(out_p, index=False)


def _change_residue_number(residue_id : str, position_change):
    """
    Residue example: "R22"
    position_change: -21
    returns: R1
    """
    aa = residue_id[0]
    num = int(residue_id[1:])
    new_num = num + position_change
    return f"{aa}{new_num}"


def _get_gipr_crop_range():
    if CROP is False:
        return slice(0, GIPR_LENGTH)
    
    return slice(21, 466-35)

def _get_glp1r_crop_range():
    if CROP is False:
        return slice(0, GLP1R_LENGTH)
    return slice(23, 463-22)



def load_reference_distances():
    ""
    script_dir = pathlib.Path(__file__).parent
    pdb_p = script_dir / 'reference/full_w_glp1_updated_chains.pdb'
    distances = get_angstrom_distances_for_pdb(pdb_p, distance_method="CA")
    distances = distances[('B', 'C')] # only keep receptors

    gipr_n_term_to_remove = int(_get_gipr_crop_range().start)
    glp1r_n_term_to_remove = int(_get_glp1r_crop_range().start)

    # crop distances
    if distances.shape[0] != GLP1R_LENGTH or distances.shape[1] != GIPR_LENGTH:
        # crop to fix issue.
        logging.warning(f'Cropping original dimer to match the cropped deep-sampling receptors.')
        logging.debug(f'Original distances: {distances.shape}')
        logging.debug(distances)
        gipr_to_keep = _get_gipr_crop_range()
        glp1r_to_keep = _get_glp1r_crop_range()
        distances = distances.iloc[glp1r_to_keep, gipr_to_keep]

        # fix residue numbering of rows and columns
        distances.index = [_change_residue_number(res_id, -glp1r_n_term_to_remove) for res_id in distances.index]
        distances.columns = [_change_residue_number(res_id, -gipr_n_term_to_remove) for res_id in distances.columns]

        logging.debug(f'Cropped distances: {distances.shape}')
        logging.debug(distances)

    # check shape again
    if distances.shape[0] != GLP1R_LENGTH or distances.shape[1] != GIPR_LENGTH:
        raise ValueError(f"Reference distances have unexpected shape: {distances.shape}. Expected (418, 410).")

    return distances


def get_contacts_intersection_percent_with_reference(model_distances : pd.DataFrame,
                                                     reference_distances : pd.DataFrame,
                                                     threshold: float = 8.0) -> int:
    """
    Returns (model contacts intersection of reference contacts) / (reference contacts) * 100
    """
    logging.debug('Model distances:')
    logging.debug(model_distances)
    logging.debug('Reference distances:')
    logging.debug(reference_distances)
    # assert cols and rows are the same
    assert model_distances.shape == reference_distances.shape, f"Distance matrices must have the same shape. Found {model_distances.shape} and {reference_distances.shape}."
    assert all(model_distances.columns == reference_distances.columns), "Distance matrices must have the same columns."
    assert all(model_distances.index == reference_distances.index), "Distance matrices must have the same index."

    # get contacts, convert df of distances to numpy arrays of bools
    model_contacts = (model_distances.values < threshold)
    reference_contacts = (reference_distances.values < threshold)
    # calculate intersection
    intersection = np.logical_and(model_contacts, reference_contacts)
    intersection_count = np.sum(intersection)
    reference_count = np.sum(reference_contacts)
    if reference_count == 0:
        return 0.0
    percent = (intersection_count / reference_count) * 100
    logging.debug(f'From {reference_count} reference contacts, found {intersection_count} shared contacts ({percent:.2f}%).')

    return percent



def add_dimer_predictions_compared_to_reference_scores(dimer_distances: pd.DataFrame,
                                                       reference_distances: pd.DataFrame) -> pd.DataFrame:
    ""
    dimer_distances = dimer_distances.copy()
    # use tqdm for progress bar, run using iterrows
    for idx, row in tqdm.tqdm(dimer_distances.iterrows(), total=dimer_distances.shape[0], desc="Comparing dimer predictions to reference"):
        dimer_dist_p = pathlib.Path(row['dimer_distances'])
        if not dimer_dist_p.exists():
            logging.warning(f"Dimer distances file does not exist: {dimer_dist_p}. Skipping.")
            continue
        model_distances = pd.read_csv(dimer_dist_p, index_col=0)
        if model_distances.shape == reference_distances.T.shape:
            logging.debug("Transposing model distances to match reference.")
            model_distances = model_distances.T
        percent_shared_contacts = get_contacts_intersection_percent_with_reference(model_distances, reference_distances)
        dimer_distances.at[idx, 'shared_contacts'] = percent_shared_contacts
    return dimer_distances


def _get_gpcr_region_dict(residue_info, slice_to_keep):
    domain_dict = dict()
    # filter for those residues in the slice_to_keep range, using 'sequence_number' column
    residue_info = residue_info.copy()
    residue_info = residue_info[residue_info['sequence_number'].between(slice_to_keep.start + 1, slice_to_keep.stop)]
    residue_info['adjusted_sequence_number'] = residue_info['sequence_number'] - slice_to_keep.start 
    for domain in residue_info['protein_segment'].unique():
        domain_residues = residue_info[residue_info['protein_segment'] == domain]
        domain_dict[domain] = tuple(domain_residues['amino_acid'] + domain_residues['adjusted_sequence_number'].astype(str))
    return domain_dict


def _get_residues_in_each_region():
    """
    Format as: 

        GLP1R: {
        'TM1' : ('N103', ..),
        'TM2' : ('..', ..),
        }
    """
    script_dir = pathlib.Path(__file__).parent
    gipr_residues = script_dir / 'domains/gipr_human_chain_info.csv'
    gipr_residues = pd.read_csv(gipr_residues)
    gipr_to_keep = _get_gipr_crop_range()
    gipr_dict = _get_gpcr_region_dict(gipr_residues, gipr_to_keep)
    
    glp1r_residues = script_dir / 'domains/glp1r_human_chain_info.csv'
    glp1r_residues = pd.read_csv(glp1r_residues)
    glp1r_to_keep = _get_glp1r_crop_range()
    glp1r_dict = _get_gpcr_region_dict(glp1r_residues, glp1r_to_keep)
    
    return gipr_dict, glp1r_dict


def add_inter_region_contacts(dimer_distances: pd.DataFrame,
                                gipr_region_dict: dict,
                                glp1r_region_dict: dict,
                                contact_threshold=8.0,
                                ):
    """
    contact_threshold: float, angstrom distance to consider a contact (between distances)
    """
    _gipr_residues = {v for values in gipr_region_dict.values() for v in values}
    _gipr_residues = sorted([sublist for sublist in _gipr_residues])
    _glp1r_residues = {v for values in glp1r_region_dict.values() for v in values}
    _glp1r_residues = sorted([sublist for sublist in _glp1r_residues])

    inter_region_pairs = itertools.product(gipr_region_dict.keys(), glp1r_region_dict.keys())
    inter_region_cols = [f'GIPR_{gipr_region}|GLP1R_{glp1r_region}|contacts' for gipr_region, glp1r_region in inter_region_pairs]
    dimer_distances.loc[:, inter_region_cols] = pd.NA

    # per dimer distance file, get the distances in subsets
    for i, row in tqdm.tqdm(dimer_distances.iterrows(), total=dimer_distances.shape[0], desc="Calculating inter-region contacts"):
        dimer_dist_p = pathlib.Path(row['dimer_distances'])
        model_distances = pd.read_csv(dimer_dist_p, index_col=0)

        for col, (gipr_region, glp1r_region) in zip(inter_region_cols, itertools.product(gipr_region_dict.keys(), glp1r_region_dict.keys())):
            gipr_residues = gipr_region_dict[gipr_region]
            glp1r_residues = glp1r_region_dict[glp1r_region]
            if not dimer_dist_p.exists():
                logging.warning(f"Dimer distances file does not exist: {dimer_dist_p}. Skipping.")
                continue
            # subset distances
            if GIPR_CHAIN == 'B' and GLP1R_CHAIN == 'C':
                sub_distances = model_distances.loc[glp1r_residues, gipr_residues]
            else:
                sub_distances = model_distances.loc[gipr_residues, glp1r_residues]
            # count contacts
            contact_count = (sub_distances.values < contact_threshold).sum()
            dimer_distances.at[i, col] = int(contact_count)


def main():
    ""
    # settings
    script_dir = pathlib.Path(__file__).parent
    log_p = script_dir / __file__.replace('.py', '.log')
    logging.basicConfig(
        filename=log_p,
        filemode='w',
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    out_p = script_dir / __file__.replace('.py', '.csv')
    residue_dist_dir = script_dir / 'residue_distances'
    residue_dist_dir.mkdir(exist_ok=True)

    # prepare distances between residues.
    logging.info("Preparing residue distances...")
    ###prepare_distances(out_p, residue_dist_dir, pdb_d=MODEL_BASE_DIR.parent / 'GLP1R_GIPR_dimers_debug_subset')
    prepare_distances(out_p, residue_dist_dir, pdb_d=MODEL_BASE_DIR, _debug=False)
    dimer_distances = pd.read_csv(out_p)

    # add inter-region distances.
    gipr_region_dict, glp1r_region_dict = _get_residues_in_each_region()
    logging.debug(f'Total residues in GIPR domains: {sum(len(v) for v in gipr_region_dict.values())}')
    logging.debug(f'Total residues in GLP1R domains: {sum(len(v) for v in glp1r_region_dict.values())}')
    add_inter_region_contacts(dimer_distances, gipr_region_dict, glp1r_region_dict)

    # load reference distances
    reference_distances = load_reference_distances()

    dimer_distances = add_dimer_predictions_compared_to_reference_scores(dimer_distances, reference_distances)
    # seperate suffix after adding the comparison scores
    dimer_distances.to_csv(out_p.with_suffix('.scores.csv'), index=False)
  

if __name__ == "__main__":
    main()
