"""

eventually:
acceleration correction when bootstrapping
"""

import os
import sys
import pathlib
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import itertools
import logging
sns.set_style("whitegrid")
sns.set_context('paper')

script_dir = pathlib.Path(__file__).parent


def glp1r(): 
    logging.info('Running GLP1R pertubation test')
    df = pd.read_csv(script_dir / 'GLP1R_jalview_annotation_fixed.csv', index_col=0)
    plot_dir = script_dir / 'pertubation_GLP1R'
    for var in VARS:
        pertubation_test(df, var, plot_dir, 'GLP1R')


def gipr():
    logging.info('Running GIPR pertubation test')
    df = pd.read_csv(script_dir / 'GIPR_jalview_annotation_fixed.csv', index_col=0)
    plot_dir = script_dir / 'pertubation_GIPR'
    
    for var in VARS:
        pertubation_test(df, var, plot_dir, 'GIPR')


def gipr_constructs():
    """
    # constructs
    # 1: 136 143 147 151
    sele gipr_construct_1, resi 136+143+147+151 and gipr
    deselect
    # 2: 138 142 146 191
    sele gipr_construct_2, resi 138+142+146+191 and gipr
    deselect
    # 3: 136 138 142 143 146 147 151 191
    sele gipr_construct_3, resi 136+138+142+143+146+147+151+191 and gipr
    deselect
    """
    df = pd.read_csv(script_dir / 'GIPR_jalview_annotation_fixed.csv', index_col=0)

    for construct_nr, construct_aas in {
        1 : [136, 143, 147, 151],
        2 : [138, 142, 146, 191],
        3 : [136, 138, 142, 143, 146, 147, 151, 191],
    }.items():
        # run pertubation test with just these residues
        plot_dir = script_dir / f'pertubation_GIPR_construct_{construct_nr}'
        construct_df = df.copy()
        # reset 'in_interface'
        construct_df['in_interface'] = False
        # only set construct_aas to in_interface = True
        construct_df['in_interface'] = construct_df['sequence_number'].isin(construct_aas)
        for var in VARS:
            pertubation_test(construct_df, var, plot_dir, f'GIPR_construct_{construct_nr}')
        construct_df.to_csv(plot_dir / f'GIPR_construct_{construct_nr}_pertubation.csv')


def glp1r_constructs():
    """
    # 1: 279 283 291
    sele glp1r_construct_1, resi 279+283+291 and glp1r
    deselect
    # 2: 283 286 290
    sele glp1r_construct_2, resi 283+286+290 and glp1r
    deselect
    # 3: 275 279 282
    sele glp1r_construct_3, resi 275+279+282 and glp1r
    deselect
    # 4: 275 279 282 283 286 290 291
    sele glp1r_construct_4, resi 275+279+282+283+286+290+291 and glp1r
    deselect
    """
    df = pd.read_csv(script_dir / 'GLP1R_jalview_annotation_fixed.csv', index_col=0)

    for construct_nr, construct_aas in {
        1 : [279, 283, 291],
        2 : [283, 286, 290],
        3 : [275, 279, 282],
        4 : [275, 279, 282, 283, 286, 290, 291],
    }.items():
        # run pertubation test with just these residues
        plot_dir = script_dir / f'pertubation_GLP1R_construct_{construct_nr}'
        construct_df = df.copy()
        # reset 'in_interface'
        construct_df['in_interface'] = False
        # only set construct_aas to in_interface = True
        construct_df['in_interface'] = construct_df['sequence_number'].isin(construct_aas)
        for var in VARS:
            pertubation_test(construct_df, var, plot_dir, f'GLP1R_construct_{construct_nr}')
        construct_df.to_csv(plot_dir / f'GLP1R_construct_{construct_nr}_pertubation.csv')


def pertubation_test(df, metric, plot_dir, gpcr):
    """ Check if the combination of residues in the predicted interface
    has higher conservation than all other pertubations.
    """
    plot_dir.mkdir(exist_ok=True)

    assert metric in df.columns, f"{metric} not in df.columns"

    tm_rows = df[df['in_TM']].index
    inter_rows = df[df['in_interface']].index
    inter_len = len(inter_rows)
    total_len = len(df)
    logging.info(f'Total residues: {total_len}')
    logging.info(f"interface residues: {inter_len}")
    logging.info(f'median {metric} of all residues: {df[metric].median()}')
    logging.info(f'median {metric} of TM residues: {df.loc[tm_rows, metric].median()}')
    logging.info(f'median {metric} of interface residues: {df.loc[inter_rows, metric].median()}')
    logging.info(f'Std {metric} of interface residues: {df.loc[inter_rows, metric].std()}')

    # plot hist with distribution of the metric in all, and of TMs
    fig, ax = plt.subplots(figsize=(4, 4))
    all_res = df[metric].values
    tm_res = df.loc[tm_rows, metric].values
    inter_res = df.loc[inter_rows, metric].values
    # inter_mean = df.loc[inter_rows, metric].mean()
    inter_median = np.median(inter_res)
    inter_sd = df.loc[inter_rows, metric].std()
    
    plt.title(f'{metric} distribution of all {gpcr} residues,\ntransmembrane residues\nand predicted dimer interface residues')
    sns.histplot(all_res, ax=ax, color='blue', label='all residues', bins=20, kde=True)
    sns.histplot(tm_res, ax=ax, color='red', label='TM residues', bins=20, kde=True)
    ax.set_xlabel(metric)
    ax.set_ylabel('count')
    plt.xlim(0, None)
    ax.legend()
    plt.savefig(plot_dir / f'distribution_{metric}_hist.png', dpi=300)
    plt.close()

    # violinplot
    fig, ax = plt.subplots(figsize=(4, 4))
    plt.title(f'{metric} distribution of all {gpcr} residues,\ntransmembrane residues\nand predicted dimer interface residues')    
    sns.violinplot(data=[all_res, tm_res, inter_res], ax=ax, color='white', alpha=0.5)
    sns.swarmplot(data=[all_res, tm_res, inter_res], ax=ax, color='black', size=2)
    ax.set_xticklabels(['all residues', 'TM residues', 'interface residues'])
    ax.set_ylabel(metric)
    plt.savefig(plot_dir / f'distribution_{metric}_violin.png', dpi=300)
    plt.close()


    # bootstrap
    num_samples = 100000
    statistics = []
    for i in range(num_samples):
        # sample random residues
        sample = np.random.choice(tm_res, inter_len, replace=True)
        statistics.append(np.median(sample))
        
    median_stat = np.median(statistics)
    std_stat = np.std(statistics)
    confidence_interval = np.percentile(statistics, [2.5, 97.5])
    logging.info(f"Median {metric} of random TM samples: {median_stat}")
    logging.info(f"Std {metric} of random TM samples: {std_stat}")
    logging.info(f"Confidence interval: {confidence_interval}")

    # plot this confidence interval and my interface values
    fig, ax = plt.subplots(figsize=(4, 4))
    plt.title(f'{metric} distribution of random {gpcr} TM samples\nand predicted dimer interface residues')
    sns.histplot(statistics, ax=ax, color='blue', label='random samples', bins=20, kde=True)
    # put the confidence interval halfway the y-axis
    conf_int_y = ax.get_ylim()[1] / 2
    # plot a confidence interval as a line with the median and 2 stds using a errorbar
    ax.errorbar(inter_median, conf_int_y, xerr=inter_sd, fmt='o', color='red', label='interface residues (median_stat, std)')
    ax.set_xlabel(metric)
    ax.set_ylabel('count')
    plt.xlim(0, None)
    ax.legend()
    plt.savefig(plot_dir / f'distribution_{metric}_bootstrap.png', dpi=300)
    plt.close()

if __name__ == "__main__":
    # set seeds
    np.random.seed(42)

    VARS = ['Conservation', 'Quality', 'Consensus', 'SHENKIN']
    # VARS= ['SHENKIN']
    script_dir = pathlib.Path(__file__).parent
    log_p = script_dir / 'step3_pertubation.log'
    logging.basicConfig(filename=log_p, level=logging.DEBUG, filemode='w')
    glp1r()
    gipr()
    glp1r_constructs()
    gipr_constructs()
