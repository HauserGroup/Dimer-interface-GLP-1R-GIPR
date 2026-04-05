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
    if plot_dir.exists():
        # if files in plot_dir, remove those
        for f in plot_dir.iterdir():
            os.remove(f)

    for var in VARS:
        pertubation_test(df, var, plot_dir, 'GLP1R', color='#F68013')


def gipr():
    logging.info('Running GIPR pertubation test')
    df = pd.read_csv(script_dir / 'GIPR_jalview_annotation_fixed.csv', index_col=0)
    plot_dir = script_dir / 'pertubation_GIPR'
    if plot_dir.exists():
        # if files in plot_dir, remove those
        for f in plot_dir.iterdir():
            os.remove(f)
    
    for var in VARS:
        pertubation_test(df, var, plot_dir, 'GIPR', color='#55A0FA')


def get_data(df, metric):
    tm_rows = df[df['in_TM']].index
    inter_rows = df[df['in_interface']].index
    inter_len = len(inter_rows)
    total_len = len(df)
    logging.info(f'Total residues: {total_len}')
    logging.info(f"interface residues: {inter_len}")
    logging.info(f'Mean {metric} of all residues: {df[metric].mean()}')
    logging.info(f'Mean {metric} of TM residues: {df.loc[tm_rows, metric].mean()}')
    logging.info(f'Mean {metric} of interface residues: {df.loc[inter_rows, metric].mean()}')
    logging.info(f'Std {metric} of interface residues: {df.loc[inter_rows, metric].std()}')
    # plot hist with distribution of the metric in all, and of TMs
    all_res = df[metric].values
    tm_res = df.loc[tm_rows, metric].values
    inter_res = df.loc[inter_rows, metric].values
    mean_inter = df.loc[inter_rows, metric].mean()
    return all_res, tm_res, inter_res, mean_inter, inter_len, tm_rows


def pertubation_test(df, metric, plot_dir, gpcr, color):
    """ Check if the combination of residues in the predicted interface
    has higher conservation than all other pertubations.
    """
    plot_dir.mkdir(exist_ok=True)

    assert metric in df.columns, f"{metric} not in df.columns"

    fig, ax = plt.subplots(figsize=(4, 4))
    all_res, tm_res, inter_res, mean_inter, inter_len, tm_rows = get_data(df, metric)
    plt.title(f'{metric} distribution of all {gpcr} residues,\ntransmembrane residues\nand predicted dimer interface residues')
    sns.histplot(all_res, ax=ax, color='blue', label='all residues', bins=20, kde=True)
    sns.histplot(tm_res, ax=ax, color='red', label='TM residues', bins=20, kde=True)
    ax.axvline(mean_inter, color='black', linestyle='--', label=f'mean interface {metric}')
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
    plt.savefig(plot_dir / f'distribution_{metric}_violin.svg')
    plt.close()


    # bootstrap
    num_samples = 100000
    statistics = []
    print(f"Sampling {num_samples} sets of size {inter_len} from {len(tm_rows)} TM residues")
    for i in range(num_samples):
        # sample random residues
        sample = np.random.choice(tm_res, inter_len, replace=True)
        statistics.append(sample.mean())
    mean_stat = np.mean(statistics)
    std_stat = np.std(statistics)
    confidence_interval = np.percentile(statistics, [2.5, 97.5])
    logging.info(f"Mean {metric} of random TM samples: {mean_stat}")
    logging.info(f"Std {metric} of random TM samples: {std_stat}")
    logging.info(f"Confidence interval: {confidence_interval}")

    # plot this confidence interval and my interface values
    sns.set_style("whitegrid")
    sns.set_context('paper')
    fig, ax = plt.subplots(figsize=(2.5, 3))
    plt.title(f'{metric} distribution of {gpcr} TM samples\nvs dimer interface\nn={num_samples} samples, size={inter_len}')
    sns.histplot(statistics, ax=ax, color=color, label='random samples', bins=20, kde=True, stat='density')
    ax.axvline(mean_inter, color='black', linestyle='--', label=f'mean interface {metric}')
    ax.axvline(confidence_interval[0], color='black', linestyle=':', label='95% CI')
    ax.axvline(confidence_interval[1], color='black', linestyle=':')
    ax.set_xlabel(metric)
    ax.set_ylabel('density')
    # set xlim to 99% of the data
    xmin = np.percentile(statistics, 0.1) - 1
    xmin = max(0, xmin)
    xmax = np.percentile(statistics, 99) + 1
    ax.set_xlim(xmin, xmax)
    ax.legend()
    plt.tight_layout()
    plt.savefig(plot_dir / f'distribution_{metric}_bootstrap.png', dpi=300)
    plt.savefig(plot_dir / f'distribution_{metric}_bootstrap.svg')
    plt.close()


def pertubation_both():
    """
    Plot bootstrap comparison of SHENKIN conservation values for interface residues
    from GLP1R and GIPR against TM-residue-based random samples.
    """
    logging.info('Running combined GLP1R and GIPR perturbation bootstrap comparison')

    df_glp1r = pd.read_csv(script_dir / 'GLP1R_jalview_annotation_fixed.csv', index_col=0)
    df_gipr = pd.read_csv(script_dir / 'GIPR_jalview_annotation_fixed.csv', index_col=0)

    metric = 'SHENKIN'
    assert metric in df_glp1r.columns and metric in df_gipr.columns, f"{metric} not found in both dataframes"

    _, tm_glp1r, inter_glp1r, mean_glp1r, len_glp1r, _ = get_data(df_glp1r, metric)
    _, tm_gipr, inter_gipr, mean_gipr, len_gipr, _ = get_data(df_gipr, metric)

    num_samples = 100000
    logging.info(f"Bootstrapping {num_samples} samples")

    stats_glp1r = [np.random.choice(tm_glp1r, len_glp1r, replace=True).mean() for _ in range(num_samples)]
    stats_gipr  = [np.random.choice(tm_gipr,  len_gipr,  replace=True).mean() for _ in range(num_samples)]

    ci_glp1r = np.percentile(stats_glp1r, [2.5, 97.5])
    ci_gipr  = np.percentile(stats_gipr,  [2.5, 97.5])

    plot_dir = script_dir / 'pertubation_combined'
    plot_dir.mkdir(exist_ok=True)

    fig, ax = plt.subplots(figsize=(3.5, 3.5))
    plt.title(f'{metric} bootstrap distribution\nGLP1R vs GIPR interface (n={num_samples})')

    sns.histplot(stats_glp1r, color='#55A0FA', label='GLP1R TM samples', bins=30, kde=True, stat='density', ax=ax)
    sns.histplot(stats_gipr, color='#F68013', label='GIPR TM samples', bins=30, kde=True, stat='density', ax=ax)

    ax.axvline(mean_glp1r, color='#55A0FA', linestyle='--', label='GLP1R interface mean')
    ax.axvline(mean_gipr, color='#F68013', linestyle='--', label='GIPR interface mean')

    ax.set_xlabel(metric)
    ax.set_ylabel('Density')
    ax.set_xlim(5, 20)
    ax.legend()
    plt.title('GIPR GLP-1R interface conservation')
    plt.tick_params(axis='both', which='both', length=0)  # remove ticks

    plt.tight_layout()
    plt.savefig(plot_dir / f'distribution_interface_{metric}_bootstrap_comparison.png', dpi=300)
    plt.savefig(plot_dir / f'distribution_interface_{metric}_bootstrap_comparison.svg')
    plt.close()


if __name__ == "__main__":
    # set seeds
    np.random.seed(42)

    # VARS = ['Conservation', 'Quality', 'Consensus', 'SHENKIN']
    VARS= ['SHENKIN']
    script_dir = pathlib.Path(__file__).parent
    log_p = script_dir / 'step3_pertubation.log'
    logging.basicConfig(filename=log_p, level=logging.DEBUG, filemode='w')
    # glp1r()
    # gipr()
    pertubation_both()
