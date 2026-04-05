import pathlib
import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns
# set fonttype none for svg to edit in illustrator
plt.rcParams['svg.fonttype'] = 'none'
sns.set_context("paper", font_scale=1.5)
sns.set_style("whitegrid")
# set default dpi to be 600 for all
plt.rcParams['figure.dpi'] = 600

COLORS = {
    'GLP1R' : '#1f77b4',       # muted blue
    'GIPR'  : '#ff7f0e',       # safety orange
    'Ambiguous' : '#7f7f7f'
}



def _plot_receptor_bound_barplot(df, save_p):
    """
    For each 'peptide' (group), plot a barplot of the fraction of models and which receptor they bind. (GLP1R, GIPR, Ambiguous)
    """
    # category
    df['peptide'] = df['peptide'].astype(str).astype('category')
    print(df['peptide'].cat.categories)
    print(df['receptor_bound'].value_counts(dropna=False))

    # barplot, where for each category, we show the fraction of models that bind to each receptor
    plt.figure(figsize=(8,4))
    # plot 3 bars side by side for each peptide
    ax = sns.histplot(
        data=df,
        x='peptide',
        hue='receptor_bound',
        multiple='dodge',
        shrink=0.8,
        palette=COLORS,
        stat='count',
        discrete=True,
        edgecolor='black',
    )

    # labels
    ax.set_xlabel("Peptide")
    # show legend with colors
    ax.set_ylabel("Number of Models")
    # add custom legend that matches colors
    handles = [plt.Line2D([0], [0], marker='o', color='w', label=key,
                          markerfacecolor=value, markersize=10) for key, value in COLORS.items()]
    ax.legend(handles=handles, title='Receptor Bound', loc='upper right')
    plt.tight_layout()
    plt.savefig(save_p)
    plt.savefig(save_p.with_suffix('.png'), dpi=600)
    plt.close()


def _plot_shared_contacts_boxplot(df, save_p):
    """
    For each 'peptide' (group), plot a boxplot of the 'shared_contacts' values.
    """
    # category
    df['peptide'] = df['peptide'].astype(str).astype('category')
    print(df['peptide'].cat.categories)

    plt.figure(figsize=(8,4))
    ax = sns.boxplot(
        data=df,
        x='peptide',
        y='shared_contacts',
        showfliers=False,
    )

    # labels
    ax.set_xlabel("Peptide")
    ax.set_ylabel("Shared Contacts")
    plt.tight_layout()
    plt.savefig(save_p)
    plt.savefig(save_p.with_suffix('.png'), dpi=600)
    plt.close()


def _plot_shared_contacts_violinplot(df, save_p):
    """
    For each 'peptide' (group), plot a violinplot of the 'shared_contacts' values.
    """
    # category
    df['peptide'] = df['peptide'].astype(str).astype('category')
    print(df['peptide'].cat.categories)

    plt.figure(figsize=(8,4))
    ax = sns.violinplot(
        data=df,
        x='peptide',
        y='shared_contacts',
        inner='box',
    )

    # labels
    ax.set_xlabel("Peptide")
    ax.set_ylabel("Shared Contacts")
    plt.tight_layout()
    plt.savefig(save_p)
    plt.savefig(save_p.with_suffix('.png'), dpi=600)
    plt.close()


def _plot_ECDF_for_shared_contacts(df, save_p):
    """
    Empirical Cumulative Distribution Function (ECDF)
    """
    # category
    df['peptide'] = df['peptide'].astype(str).astype('category')
    print(df['peptide'].cat.categories)

    plt.figure(figsize=(8,4))
    ax = sns.ecdfplot(
        data=df,
        x='shared_contacts',
        hue='peptide',
        palette='Set2',
    )

    # labels
    ax.set_xlabel("Shared Contacts")
    ax.set_ylabel("ECDF")
    plt.tight_layout()
    plt.savefig(save_p)
    plt.savefig(save_p.with_suffix('.png'), dpi=600)
    plt.close()


def _stripplot_shared_contacts(df, save_p):
    """
    Seaborn swarmplot
    """
    # category
    df['peptide'] = df['peptide'].astype(str).astype('category')

    plt.figure(figsize=(8,4))
    ax = sns.stripplot(
        data=df,
        y='shared_contacts',
        x='peptide',
        palette='Set2',
        size=5,        
    )

    # labels
    ax.set_ylabel("Shared Contacts")
    ax.set_xlabel("peptide")
    ax.set_ylim(0, 100)

    plt.tight_layout()
    plt.savefig(save_p)
    plt.savefig(save_p.with_suffix('.png'), dpi=600)
    plt.close()


def _plot_contact_per_region(df, peptide, save_p, only_use_TMs=False):

    order=['N-term', 'TM1', 'ICL1', 'TM2', 'ECL1', 'TM3', 'ICL2', 'TM4', 'ECL2', 'TM5', 'ICL3', 'TM6', 'ECL3', 'TM7', 'H8', 'C-term']
    
    df = df.copy()
    df = df[df['peptide'] == peptide]

    contact_cols = [c for c in df.columns if c.endswith('|contacts')]
    gipr_regions = [c.split('|')[0] for c in contact_cols]
    gipr_regions = [c.replace('GIPR_', '') for c in gipr_regions]

    glp1r_regions = [c.split('|')[1] for c in contact_cols]
    glp1r_regions = [c.replace('GLP1R_', '') for c in glp1r_regions]

    # fill a heatmap
    plt.figure(figsize=(10,6))
    heatmap_data = []
    for gipr_r, glp1r_r, contact_col in zip(gipr_regions, glp1r_regions, contact_cols):
        avg_contact = df[contact_col].mean()
        heatmap_data.append({
            'GIPR_Region': gipr_r,
            'GLP1R_Region': glp1r_r,
            'Avg_Contacts': avg_contact
        })

    heatmap_df = pd.DataFrame(heatmap_data)
    heatmap_pivot = heatmap_df.pivot(index='GIPR_Region', columns='GLP1R_Region', values='Avg_Contacts')
    # sort index and columns on order
    heatmap_pivot = heatmap_pivot.reindex(index=order, columns=order)
    if only_use_TMs:
        tm_order = [r for r in order if r.startswith('TM')]
        heatmap_pivot = heatmap_pivot.reindex(index=tm_order, columns=tm_order)

    # <0.1 mask
    mask = heatmap_pivot < 1.0
    ax = sns.heatmap(
        heatmap_pivot,
        annot=True,
        fmt=".2f",
        cmap="YlGnBu",
        cbar_kws={'label': 'Average Contacts'},
        mask=mask,
    )
    # title
    ax.set_title(f"Average Contacts per Region for {peptide}\nMasked < 1.0")
    ax.set_xlabel("GLP1R Region")
    ax.set_ylabel("GIPR Region")
    plt.tight_layout()
    plt.savefig(save_p)
    plt.savefig(save_p.with_suffix('.png'), dpi=600)
    plt.close()


def main():
    """
    stem,seed,peptide,model,sample,iptm,id,pdb_p,dimer_distances,GLP1R_pep_distances,GIPR_pep_distances,receptor_bound,shared_contacts
    """
    script_dir = pathlib.Path(__file__).parent.resolve()
    scores = pd.read_csv(script_dir / "run.scores.csv")

    for (plot_dir, take_best) in zip(
        [script_dir / 'plots_all_models',  script_dir / 'plots_best_per_peptide_seed_model'],
        [False, True]
    ):
        df = scores.copy()
        plot_dir.mkdir(exist_ok=True)
        take_best = True
        if take_best:
            # keep best iptm per (peptide, seed)
            df = df.sort_values('iptm', ascending=False).groupby(['peptide', 'seed', 'model']).head(1).reset_index(drop=True)

        for pep in df['peptide'].unique():
            _plot_contact_per_region(df.copy(), peptide=pep, save_p=plot_dir / f"contact_per_region_{pep}.svg", only_use_TMs=True)
            _plot_contact_per_region(df.copy(), peptide=pep, save_p=plot_dir / f"contact_per_region_{pep}_all_regions.svg", only_use_TMs=False)
        _plot_receptor_bound_barplot(df.copy(), plot_dir / "receptor_bound_barplot.svg")
        _plot_shared_contacts_boxplot(df.copy(), plot_dir / "shared_contacts_boxplot.svg")
        _plot_shared_contacts_violinplot(df.copy(), plot_dir / "shared_contacts_violinplot.svg")
        _plot_ECDF_for_shared_contacts(df.copy(), plot_dir / "shared_contacts_ECDF.svg")
        _stripplot_shared_contacts(df.copy(), plot_dir / "shared_contacts_stripplot.svg")


if __name__ == "__main__":
    main()
