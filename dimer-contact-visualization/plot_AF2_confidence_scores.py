import pandas as pd
import pathlib
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap


def af2_cmap():
    """ These are the descriptions of the iptm+ptm values.

    100 to 90 – high accuracy expected
    90 to 70 – backbone expected to be modeled well
    70 to 50 – low confidence, caution
    50 to 0 – should not be interpreted, may be disordered 
    """
    # create the range groups
    range_groups = [(0, 50), (50, 70), (70, 90), (90, 100)]
    # bad = red to orange
    bad = [(1, 0, 0), (1, 0.5, 0)]
    # low = orange to yellow
    low = [(1, 0.5, 0), (1, 1, 0)]
    # medium = yellow to light blue
    medium = [(1, 1, 0), (0, 0, 1)]
    # high = light blue to dark blue
    high = [(0, 0, 1), (0, 0, 0.5)]

    # create the color map
    colors = bad + low + medium + high
    segments = 200  # Number of color segments
    cmap_name = 'AlphaFold2_ColorMap'
    cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=segments,
                                           )
                                           
    return cm


def plot_heatmap_of_scores(df, save_path, score_col="iptm_ptm_score_of_best"):
    """
    use seaborn to visualize the scores with a heatmap.
    # the top row will be the peptides, with:
        "No peptide", "gip", "glp1"
    # the left column will be the dimer pairs:
        "glp1r x glp1r", "glp1r x gipr", "gipr x gipr"

    # the heatmap will be the iptm+ptm score
    """
    # heatmap df
    column_headers = ["no_peptide", "gip", "glp1"]
    row_headers = ["glp1r x glp1r", "glp1r x gipr", "gipr x gipr", ]
    # create an empty df with the correct columns and rows
    df_pivot = pd.DataFrame(columns=column_headers, index=row_headers)
    # fill in the df with the iptm+ptm scores
    for index, row in df.iterrows():
        if row["peptide"] is np.NaN:
            peptide = "no_peptide"
        else:
            peptide = row["peptide"]
        column = row["gpcr_a"].replace("_human", "") + " x " + row["gpcr_b"].replace("_human", "")
        df_pivot.at[column, peptide] = row[score_col]
    df_pivot = df_pivot.astype(float)

    # set the max and min values
    vmin = 0
    vmax = 1
    sns.heatmap(df_pivot, annot=True, 
                # use the alphafold2 color map
                cmap=af2_cmap(),
                # put the column names at the top
                cbar_kws={"orientation": "horizontal",
                           "label": "IPTM+PTM score",
                           "shrink": 0.5,
                          "ticks": [0, 0.25, 0.5, 0.75, 1],
                            },
                # set the min and max values
                vmin=vmin, 
                vmax=vmax,
                # set the font size
                annot_kws={"size": 10,
                           })
    
    plt.xlabel("Peptide")
    plt.ylabel("GPCR dimer")
    # rotate the y labels
    plt.yticks(rotation=0)
    plt.title(f"{score_col} scores of GPCR dimers")
    # make it fit
    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()

if __name__ == "__main__":
    # plot pLDDT
    PLDDT_DF = pd.read_csv(pathlib.Path("data/dimer_pLDDT_scores.csv"))
    plot_heatmap_of_scores(PLDDT_DF,
                           "plots/score_heatmaps/pLDDT_heatmap.png", 
                           score_col="pLDDT")
    
    # plot iptm+ptm
    IPTM_PTM_DF = pd.read_csv(pathlib.Path("data/dimer_iptm_ptm_scores.csv"))
    plot_heatmap_of_scores(IPTM_PTM_DF,
                           "plots/score_heatmaps/iptm_ptm_heatmap.png", 
                           score_col="iptm_ptm_score_of_best")
    
