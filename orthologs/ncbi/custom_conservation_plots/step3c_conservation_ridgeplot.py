"""
Get residues positions for:
- GLP1R interface
- GLP1R clusters
- GIPR interface
- GIPR clusters

For these residue groups, plot distributions of SHENKIN in ridgeplots
"""
import pandas as pd
import pathlib
from step3_pertubation import get_data


import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns
import pathlib
import sys 
import numpy as np
sys.path.append(str(pathlib.Path(__file__).parent.parent))
from joypy import joyplot
import ast
from matplotlib.colors import ListedColormap, to_hex

# style
sns.set_context("paper")
sns.set_style("whitegrid")
sns.set_style("ticks", {"xtick.major.pad": 0, "ytick.major.pad": 0})
plt.rcParams.update({
    "font.size": 8,
    "axes.titlesize": 10,
    "axes.labelsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "figure.titlesize": 10,
})



def plot_ridgeplot_from_df(
    df: pd.DataFrame,
    save_path: pathlib.Path,
    x_range: tuple = (0.0, 1.0),
    y_col : str = 'boot_mean_identity',
    overlap: float = 0.65,
    ylim: float = None,
):
    """
    Horizontal ridgeline (joy) plot of boot_mean_identity by receptor family.

    Parameters
    ----------
    df : DataFrame
        Needs columns: 'boot_mean_identity', 'name', 'color'
    save_path : pathlib.Path
        PNG is written here; an SVG with the same stem is written too.
    x_range : tuple, default (0.0, 1.0)
        Limits for the x-axis (identity runs 0→1).
    y_col : str = 'boot_mean_identity',
        Column name in `df` to plot on the y-axis.
    overlap : float, default 0.65
        0 = no overlap, 1 = full overlap.
    """
    # ---------------------------------------------------------------- set-up
    df = df.copy()

    # set colors
    groups = df['name'].unique().tolist()
    name_to_col = df.set_index('name')['color'].to_dict()
    colors = [name_to_col[g] for g in groups]
    colors = [to_hex(c) for c in colors]
    # sort df on name by groups
    df['name'] = df['name'].astype('category')
    df['name'] = df['name'].cat.reorder_categories(groups, ordered=True)
    print(f"Colors: {colors}")
    print(f"Groups: {groups}")
    # sort so it is 
    # ---------------------------------------------------------------- plot
    fig, axes = joyplot(
        data=df,
        by="name",
        column=y_col,
        #fade=True,
        ylim=ylim,
        range_style='all', # own, group, all
        grid=True,
        overlap=overlap,
        figsize=(2, 2),
        x_range=x_range,
        linewidth=0.25,
        linecolor='black',
        alpha=1.0,
        ylabelsize=4,
        color=colors,
    )

    # plot xticks
    axes[-1].set_xticks(np.arange(x_range[0], x_range[1] + 0.01, 10))
    axes[-1].tick_params(axis='x', labelsize=6)
        
    # Cosmetics: axis labels & spines
    axes[-1].set_xlabel(y_col)
    for spine in ("left", "right", "top"):
        axes[-1].spines[spine].set_visible(False)

    # ------------------------------------------------------------- export
    fig.savefig(save_path, dpi=300)
    fig.savefig(save_path.with_suffix(".svg"))
    plt.close(fig)


# settings
script_dir = pathlib.Path(__file__).parent
metric = 'SHENKIN'
df_glp1r = pd.read_csv(script_dir / 'GLP1R_jalview_annotation_fixed.csv', index_col=0)
df_gipr = pd.read_csv(script_dir / 'GIPR_jalview_annotation_fixed.csv', index_col=0)

glp1r_clusters = {
    'cluster 6': [223, 226, 227, 235, 236, 237, 240, 243, 244, 245, 247, 248, 250, 251, 254, 255],
    'cluster 7': [187, 188, 190, 191, 205],
    'cluster 8': [179, 181]
}
glp1r_clusters_colors = {
    'cluster 6': 'lightcyan',
    'cluster 7': 'slateblue',
    'cluster 8': 'purple',
}
gipr_clusters = {
    'cluster 1': [265, 267],
    'cluster 2': [144, 145, 148, 155, 163],
    'cluster 3': [84, 88],
    'cluster 4': [68, 70, 71],
}
gipr_clusters_colors = {
    'cluster 1': 'darkkhaki',
    'cluster 2': 'yellow',
    'cluster 3': 'black',
    'cluster 4': 'orange',
}

# report for whole interface
glp1r_tms = df_glp1r[df_glp1r['in_TM']][metric].values
glp1r_interface = df_glp1r[df_glp1r['in_interface']][metric].values
print(f'Mean glp1r_interface {metric}: {glp1r_interface.mean()}')
gipr_tms = df_gipr[df_gipr['in_TM']][metric].values
gipr_interface = df_gipr[df_gipr['in_interface']][metric].values
print(f'Mean gipr_interface {metric}: {gipr_interface.mean()}')


# report for clusters
for clust, aas in glp1r_clusters.items():
    cluster_data = df_glp1r[df_glp1r['sequence_number'].isin(aas)][metric].values
    print(f'Mean GLP1R {clust} {metric}: {cluster_data.mean()}')
    print(f"Values for GLP1R {clust} {metric}: {', '.join(map(str, cluster_data))}")
for clust, aas in gipr_clusters.items():
    cluster_data = df_gipr[df_gipr['sequence_number'].isin(aas)][metric].values
    print(f'Mean GIPR {clust} {metric}: {cluster_data.mean()}')
    print(f"Values for GIPR {clust} {metric}: {', '.join(map(str, cluster_data))}")

# merge all into a dataframe
data = []
data.append(('GIPR TMs', gipr_tms, 'peachpuff'))
data.append(('GIPR interface',  gipr_interface, 'darkorange'))
for clust, aas in gipr_clusters.items():
    cluster_data = df_gipr[df_gipr['sequence_number'].isin(aas)][metric].values
    cluster_color = gipr_clusters_colors[clust]
    data.append((f'GIPR {clust}', cluster_data, cluster_color))


data.append(('GLP1R TMs',  glp1r_tms, 'lightskyblue'))
data.append(('GLP1R interface',  glp1r_interface, 'royalblue'))
for clust, aas in glp1r_clusters.items():
    cluster_data = df_glp1r[df_glp1r['sequence_number'].isin(aas)][metric].values
    cluster_color = glp1r_clusters_colors[clust]
    data.append((f'GLP1R {clust}', cluster_data, cluster_color))


# create dataframe
df = pd.DataFrame(data, columns=['name', metric, 'color'])
# explode the values column
df = df.explode(metric)
# convert values to float
df[metric] = df[metric].astype(float)
df.to_csv(script_dir / 'step3c_ridgeplot_shenkin.csv', index=False)
# duplicate the metric lists
df = df.copy()

# plot
save_path = script_dir / 'step3c_ridgeplot_shenking.png'
# reverse the df order (by name) from the original df
current_order = df['name'].unique().tolist()
# reverse the order
df['name'] = pd.Categorical(df['name'], categories=current_order[::-1], ordered=True)
# cluster 3 is not visible (only 2 vals)
# add 1 more datapoint selected so it appears as a line
df = df.copy()
# insert row for GIPR cluster 3
gipr_clust3_mean = df[df['name'] == 'GIPR cluster 3'][metric].mean()
gipr_clust3_color = gipr_clusters_colors['cluster 3']
df = pd.concat([
    df,
    pd.DataFrame(
        {'name': ['GIPR cluster 3'], f'{metric}' : [gipr_clust3_mean+0.1], 'color': gipr_clust3_color},
        {'name': 'GIPR cluster 3', metric: gipr_clust3_mean-0.1, 'color': gipr_clust3_color}
    )
], ignore_index=True)
plot_ridgeplot_from_df(
    df,
    save_path=save_path,
    x_range=(0.0, 30.0),
    y_col=metric,
    overlap=0.0,
)