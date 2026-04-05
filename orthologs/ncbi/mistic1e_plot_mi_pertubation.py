import pandas as pd 
import pathlib
import os
import matplotlib.pyplot as plt
import seaborn as sns
import random
import numpy as np



def run_pertubation_test(df, plot_p=False, samples=10_000, sample_size=None, random_state=None):
    """
    Compare mean MI values of interface vs. non-interface residues using a one-sided
    permutation test (H₀: means are equal; H₁: interface > non-interface).

    Parameters
    ----------
    df : pandas.DataFrame
        Must contain columns:
        - 'MIvalue'   : numeric mutual-information scores
        - 'interface' : any non-null marker for interface residues
    plot_p : str or bool, default False
        If a string (path), save a histogram of the permutation distribution with
        the observed statistic marked to this file. If True, show the plot.
    samples : int, default 10_000
        Number of random permutations.
    sample_size : int or None, default None
        Number of residues to include from *each* group when computing the
        statistic. If None (recommended), the full group sizes are used. If given,
        both groups are down-sampled without replacement to this size on every
        iteration.
    random_state : int or None
        Seed for reproducibility.

    Returns
    -------
    p_value : float
        One-sided p-value (fraction of permuted differences ≥ observed).
    observed_diff : float
        Observed difference in group means (interface − non-interface).
    perm_diffs : ndarray
        Array of length `samples` with the permuted mean differences.

    Notes
    -----
    * Applies a +1 / (n+1) correction to avoid zero p-values.
    * If `plot_p` is a filepath, the plot is saved; if True, it is shown interactively.
    """
    rng = np.random.default_rng(random_state)

    # --- Pull MI values from labeled groups ---
    int_vals  = df.loc[df['interface'].notna(), 'MIvalue'].dropna().to_numpy()
    non_vals  = df.loc[df['interface'].isna(),  'MIvalue'].dropna().to_numpy()

    if int_vals.size == 0 or non_vals.size == 0:
        raise ValueError("Both groups must contain at least one non-NaN MI value.")

    # --- Compute observed statistic ---
    observed_diff = int_vals.mean() - non_vals.mean()

    # --- Prepare for permutation ---
    group_sizes = (int_vals.size, non_vals.size)
    pooled = np.concatenate([int_vals, non_vals])
    perm_diffs = np.empty(samples)

    downsample = (
        sample_size is not None and
        sample_size <= min(group_sizes)
    )

    # --- Permutation loop ---
    for i in range(samples):
        rng.shuffle(pooled)
        grp_a = pooled[:group_sizes[0]]
        grp_b = pooled[group_sizes[0]:]

        if downsample:
            grp_a = rng.choice(grp_a, sample_size, replace=False)
            grp_b = rng.choice(grp_b, sample_size, replace=False)

        perm_diffs[i] = grp_a.mean() - grp_b.mean()

    # --- Compute one-sided p-value ---
    p_value = (np.sum(perm_diffs >= observed_diff) + 1) / (samples + 1)

    # --- Optional plot ---
    if plot_p:
        plt.figure(figsize=(8, 5))
        plt.hist(perm_diffs, bins=50, alpha=0.75, edgecolor='black')
        plt.axvline(observed_diff, color='red', linestyle='--', lw=2,
                    label=f"Observed = {observed_diff:.4f}")
        plt.title(f"Permutation test\nOne-sided p = {p_value:.4g}")
        plt.xlabel("Δ mean MI (interface − non-interface)")
        plt.ylabel("Frequency")
        plt.legend()
        plt.tight_layout()

        if isinstance(plot_p, str) or isinstance(plot_p, pathlib.Path):
            plt.savefig(plot_p)
            plt.close()
        else:
            plt.show()

    return p_value, observed_diff, perm_diffs


def get_glp1r_interface_residues():
    return [279,289,282,283,286,287,290,291]


def get_gipr_interface_residues():
    return [136,147,192,150,151,146,142,143,188,139,191,138,141]

# paths
script_dir = pathlib.Path(__file__).parent
out_dir = script_dir / 'OLD_MISTIC1_MI_pertubation'
out_dir.mkdir(exist_ok=True)
df_p = out_dir / 'MI_annotated_interface.csv'

print("loading or creating dataframe...")
# define which pairs belong to an interface
if not os.path.exists(df_p):
    glp1r_interface_residues = get_glp1r_interface_residues()
    gipr_interface_residues = get_gipr_interface_residues()
    df = pd.read_csv(script_dir / 'MISTIC1_MI/MI_annotated.csv')
    df['interface'] = pd.NA
    # match gpcr_sequence_number and 'GPCR1' to set interface status (pd.NA or 'GLP1R' or 'GIPR')
    for i, row in df.iterrows():
        if row['GPCR1'] == 'GLP1R':
            if row['gpcr_sequence_number1'] in glp1r_interface_residues:
                df.at[i, 'interface'] = 'GLP1R'
        elif row['GPCR1'] == 'GIPR':
            if row['gpcr_sequence_number1'] in gipr_interface_residues:
                df.at[i, 'interface'] = 'GIPR'

    # match gpcr_sequence_number and 'GPCR2' to set interface status (pd.NA or 'GLP1R' or 'GIPR')
    for i, row in df.iterrows():
        old_interface = row['interface']
        if row['GPCR2'] == 'GLP1R':
            if row['gpcr_sequence_number2'] in glp1r_interface_residues:
                if pd.isna(old_interface):
                    df.at[i, 'interface'] = 'GLP1R'
                else:
                    df.at[i, 'interface'] = f'{old_interface}|GLP1R'
        elif row['GPCR2'] == 'GIPR':
            if row['gpcr_sequence_number2'] in gipr_interface_residues:
                df.at[i, 'interface'] = 'GIPR'
                if pd.isna(old_interface):
                    df.at[i, 'interface'] = 'GIPR'
                else:
                    df.at[i, 'interface'] = f'{old_interface}|GIPR'

    df.to_csv(df_p, index=False)
else:
    print('df already exists')
    df = pd.read_csv(df_p)
if 'interface' not in df.columns:
    raise ValueError('Interface column not found in dataframe')
print(df)


# # plot overall conservation pertubation
# for (sample_size, samples) in [(1, 10000), (5, 10000), (10, 10000)]:
#     out = run_pertubation_test(df, 
#                          out_dir / f'pertubation_test_{sample_size}size_{samples}samples.png', 
#                          samples=samples, 
#                          sample_size=sample_size)
#     print(f"Sample size: {sample_size}, Samples: {samples},\n{out}")

# filter for only TM-TM interactions
df['IntraPair'] = df['GPCR1'] == df['GPCR2']
df = df[df['segment1'].str.upper().str.contains('TM') & df['segment2'].str.upper().str.contains('TM')]

# plot IntraPair only
intrapair = df[df['IntraPair']]
# both gpcrs are same
intrapair_glp1r = intrapair[intrapair['GPCR1'] == 'GLP1R']
intrapair_gipr = intrapair[intrapair['GPCR1'] == 'GIPR']
# plot permutation for intrapair (5, 10000)
out = run_pertubation_test(intrapair_glp1r,
                            out_dir / f'intrapair_pertubation_test_glp1r.png',
                            samples=10000,
                            sample_size=5)
print(f"GLP1R: Sample size: {5}, Samples: {10000},\n{out}")
out = run_pertubation_test(intrapair_gipr,
                            out_dir / f'intrapair_pertubation_test_gipr.png',
                            samples=10000,
                            sample_size=5)
print(f"GIPR: Sample size: {5}, Samples: {10000},\n{out}")


# keep only InterPairs
df = df[~df['IntraPair']] 
df.to_csv(out_dir / 'MI_inter_annotated_TM.csv', index=False)


# for (sample_size, samples) in [(1, 10000), (5, 10000), (10, 10000)]:
#     out = run_pertubation_test(df, 
#                          out_dir / f'inter_pertubation_test_{sample_size}size_{samples}samples.png', 
#                          samples=samples, 
#                          sample_size=sample_size)
#     print(f"Sample size: {sample_size}, Samples: {samples},\n{out}")


# plot violinplot distributions of non-interface GIPR vs GLP1R vs GLP1R|GIPR
subplots = 4
fig, axes = plt.subplots(ncols=subplots, figsize=(8, 4), sharey=True)
plt.subplots_adjust(hspace=0.5)
df['interface'] = df['interface'].fillna('')
# plot violinplot of all pairs
sns.violinplot(y='MIvalue', data=df[df['interface'] == ''], ax=axes[0], inner='quartile')
axes[0].set_title('Non-interface')
sns.violinplot(y='MIvalue', data=df[df['interface'].str.contains('GLP1R')],
               ax=axes[1], inner='quartile')
axes[1].set_title('GLP1R')
sns.violinplot(y='MIvalue', data=df[df['interface'].str.contains('GIPR')],
               ax=axes[2], inner='quartile')
axes[2].set_title('GIPR')
sns.violinplot(y='MIvalue', data=df[df['interface'] == 'GLP1R|GIPR'], ax=axes[3], inner='quartile')
axes[3].set_title('GLP1R|GIPR')
plt.tight_layout()
plt.savefig(out_dir / 'MI_annotated_interface_violinplot.png', dpi=600)
plt.close()