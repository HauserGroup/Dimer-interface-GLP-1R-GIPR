""" 
"""
import pathlib
import pandas as pd
import seaborn as sns 
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter


def get_tm_only_df():
    script_dir = pathlib.Path(__file__).parent
    out_dir = script_dir / 'MISTIC1_MI_pertubation'
    df = pd.read_csv(out_dir / 'MI_inter_annotated_constructs.csv')

    # only keep pairs involving two TMs
    # TM in segment1.upper() AND segment2.upper()
    print("Filtering pairs with two TM segments")
    print(f"Original dataframe shape: {df.shape}")
    df = df[df['segment1'].str.upper().str.contains('TM') & df['segment2'].str.upper().str.contains('TM')]
    print(f"Filtered dataframe shape: {df.shape}")
    return df


def construct_vs_non_construct_MI():
    """
    Plot conservation of regular TM inter-receptor pairs.
    Compared to pairs considered for the interface (Construct != None)
    """
    df = get_tm_only_df()
    df['IntraPair'] = df['GPCR1'] == df['GPCR2']
    # keep only InterPairs
    df = df[~df['IntraPair']] 

    script_dir = pathlib.Path(__file__).parent
    out_dir = script_dir / 'MISTIC1_MI_pertubation'

    number_of_constructs = list(range(1, 8))

    # part of construct if 'Construct' not nan
    df['PartOfConstruct'] = df['Construct'].notna().astype(int)

    # plot violin plot conservation
    plt.figure(figsize=(4, 6))
    # calculate difference in median
    median = df.groupby('PartOfConstruct')['MIvalue'].median()
    diff = median[1] - median[0]

    # plot violin
    sns.violinplot(y='MIvalue', hue='PartOfConstruct', data=df, split=True, inner='quartile')
    difference_in_median_label = f'Difference in median\n(Non-costruct - construct)\n= {diff:.2f}'
    plt.title(difference_in_median_label, fontsize=10)
    plt.suptitle('Conservation of inter GIPR-GLP1R TM pairs',
                 # center
                    x=0.5)
    plt.tight_layout()
    plt.savefig(out_dir / 'TM_pairs_conservation.png', dpi=600)
    plt.close()
    


def non_construct_vs_individual_constructs_MI():
    """
    """
    df = get_tm_only_df()
    df['PartOfConstruct'] = df['Construct'].notna().astype(int)
    df['Construct'] = df['Construct'].fillna('')
    # check if GPCR1 and GPCR2 are the same
    df['IntraPair'] = df['GPCR1'] == df['GPCR2']
    # keep only InterPairs
    df = df[~df['IntraPair']] 

    script_dir = pathlib.Path(__file__).parent
    out_dir = script_dir / 'MISTIC1_MI_pertubation'
    number_of_constructs = list(range(1, 8))
    number_of_violinplots = len(number_of_constructs) + 1

    fig, axes = plt.subplots(ncols=number_of_violinplots, figsize=(8, 3), sharey=True)
    fig.subplots_adjust(hspace=0.5)
    fig.suptitle('Conservation of inter GIPR-GLP1R TM pairs', fontsize=10)

    # example, 'Construct' == 6|7, therefore, part of construct 6 and 7
    for construct_number in number_of_constructs:
        construct_df = df[df['Construct'].str.contains(f'{construct_number}')]
        sns.violinplot(y='MIvalue', hue='PartOfConstruct', data=construct_df,
                       split=False, inner='quartile', ax=axes[construct_number],
                       legend=False)
        axes[construct_number].set_title(f'M{construct_number}')

    # plot non-construct pairs
    sns.violinplot(y='MIvalue', hue='PartOfConstruct', data=df[df['PartOfConstruct'] == 0],
                   split=False, inner='quartile', ax=axes[0],
                   legend=False)
    axes[0].set_title('Non-construct')
    # set ymax to 10
    axes[0].set_ylim(-10, 10)
    plt.tight_layout()
    plt.savefig(out_dir / 'TM_pairs_conservation_individual_constructs.png', dpi=600)
    plt.close()


def plot_pertubation():
    """
    """
    n_samples = 10000
    sample_size = 5
    out_dir = pathlib.Path(__file__).parent / 'MISTIC1_MI_pertubation'
    out_dir.mkdir(exist_ok=True)
    random_seed = 42
    rng = np.random.default_rng(random_seed)

    df = get_tm_only_df()
    df['IntraPair'] = df['GPCR1'] == df['GPCR2']
    # keep only InterPairs
    df = df[~df['IntraPair']]
    # sort on evalue, high to low
    df = df.sort_values(by='MIvalue', ascending=False)
    df.to_csv(out_dir / 'only_inter_tm_df.csv', index=False)

    # Observed difference in MIvalue between with and without Construct
    has_construct = df[~df['Construct'].isna()]
    no_construct = df[df['Construct'].isna()]
    observed_diff = has_construct['MIvalue'].mean() - no_construct['MIvalue'].mean()

    # Permutation test
    has_construct = df[~df['Construct'].isna()]['MIvalue'].dropna()
    no_construct  = df[df['Construct'].isna()]['MIvalue'].dropna()
    if len(has_construct) < sample_size or len(no_construct) < sample_size:
        raise ValueError("Not enough samples in one of the groups for the requested sample size.")

    # Observed difference based on full data
    observed_diff = has_construct.mean() - no_construct.mean()

    # Permutation test: subsample and compute difference
    null_diffs = np.empty(n_samples)
    for i in range(n_samples):
        sample1 = rng.choice(has_construct, size=sample_size, replace=False)
        sample2 = rng.choice(no_construct, size=sample_size, replace=False)
        null_diffs[i] = sample1.mean() - sample2.mean()

    # Two-sided p-value
    p_value = (np.sum(np.abs(null_diffs) >= np.abs(observed_diff)) + 1) / (n_samples + 1)

    # Save outputs
    np.save(out_dir / 'null_distribution.npy', null_diffs)
    with open(out_dir / 'observed.txt', 'w') as f:
        f.write(f"Observed difference (full data): {observed_diff:.4f}\n")
        f.write(f"Permutation P-value (n={n_samples}, subsample size={sample_size}): {p_value:.4f}\n")

    # -------------------------------------------------------------------------
    # Plotting with scientific style
    # -------------------------------------------------------------------------
    plt.figure(figsize=(5.5, 5))
    plt.hist(null_diffs, bins=25, color='steelblue', alpha=0.75, edgecolor='black')
    plt.axvline(observed_diff, color='crimson', linestyle='--', linewidth=2)

    # Add text annotation inside the plot
    legend_labels = [
        'Empirical null\ndistribution',
        f'Observed Δ =\n{observed_diff:.4f}',
    ]
    # also add p-value
    legend_labels.append(f'Permutation P-value = {p_value:.4f}')
    # visualize p-value
    plt.text(
        x=observed_diff + 0.5,
        y=plt.ylim()[1] * 0.96,
        s=f'p = {p_value:.4f}',
        fontsize=10,
        color='black'
    )

    plt.legend(legend_labels, fontsize=10, loc='upper right', frameon=False)

    # Axis labels and style
    plt.xlabel('Δ mean MIvalue (sample of 5)', fontsize=12)
    plt.ylabel('Frequency', fontsize=12)
    plt.title(
        f'Permutation Test on MIvalue\n(Observed vs Null Distribution)',
        fontsize=13
    )

    # Format y-axis with scientific notation
    ax = plt.gca()
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 3))
    ax.tick_params(labelsize=10)
    plt.tight_layout()
    plt.savefig(out_dir / 'permutation_test_plot.png', dpi=300)
    plt.savefig(out_dir / 'permutation_test_plot.svg')
    plt.close()
    print(f"Permutation test complete. Results saved to {out_dir.resolve()}")


if __name__ == "__main__":
    # construct_vs_non_construct_MI()
    # non_construct_vs_individual_constructs_MI()
    plot_pertubation()

