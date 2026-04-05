"""
1. match proteins by species.

"""
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pathlib
from Bio import SeqIO


def get_tm_cols(cols, tm_residues):
    tm_cols = []
    for col in cols:
        res_part = col.partition('_')[2]
        seq_number = res_part[1:]
        assert seq_number.isdigit()
        if int(seq_number) in tm_residues:
            tm_cols.append(col)
    return tm_cols


def get_gipr_residue_to_chain():
    df = pd.read_csv('code/gipr_human_chain_info.csv')
    #     sequence_number,amino_acid,protein_segment,display_generic_number
    # 1,M,N-term,
    # 2,T,N-term,
    df['residue'] = 'gipr_' + df['amino_acid'] + df['sequence_number'].astype(str)
    return df.set_index('residue')['protein_segment'].to_dict()


def get_glp1r_residue_to_chain():
    df = pd.read_csv('code/glp1r_human_chain_info.csv')
    df['residue'] = 'glp1r_' + df['amino_acid'] + df['sequence_number'].astype(str)
    return df.set_index('residue')['protein_segment'].to_dict()


if __name__ == "__main__":
    script_dir = pathlib.Path(__file__).parent.absolute()
    df = pd.read_csv(script_dir / 'step5_glp1r_gipr_co_mutations.csv', index_col=0)
    glp1r_data = pd.read_csv(script_dir / 'GLP1R_jalview_annotation_fixed.csv', index_col=0)
    glp1r_tm_residues = list(set(glp1r_data[glp1r_data['in_TM']]['sequence_number'].values))
    glp1r_tm_cols = get_tm_cols(df.index, glp1r_tm_residues)
    glp1r_residue_to_chain = get_glp1r_residue_to_chain()
    glp1r_tm_index_groups = [glp1r_residue_to_chain[i] for i in glp1r_tm_cols]
    print(glp1r_tm_index_groups)
    gipr_data = pd.read_csv(script_dir / 'GIPR_jalview_annotation_fixed.csv', index_col=0)
    gipr_tm_residues = gipr_data[gipr_data['in_TM']]['sequence_number'].values
    gipr_tm_cols = list(set(get_tm_cols(df.columns, gipr_tm_residues)))
    gipr_residue_to_chain = get_gipr_residue_to_chain()
    gipr_tm_index_groups = [gipr_residue_to_chain[i] for i in gipr_tm_cols]
    
    df_only_tm = df.loc[glp1r_tm_cols, gipr_tm_cols]
    df_only_tm.to_csv(script_dir / 'step5b_glp1r_gipr_co_mutations_tm.csv')

    # create seaborn heatmap
    sns.set(font_scale=0.5)
    sns.set_style("whitegrid")
    # Create a heatmap
    sns.heatmap(
        df,
        cmap="coolwarm",  # Color map for better visualization
        annot=False,  # Show cell values
        # fmt=".2f",  # Format cell values
        cbar=True,  # Show color bar
        square=True,  # Make cells square
        # deal with nans
        mask=np.isnan(df),
        # color by group
        row_colors=glp1r_tm_index_groups,
        col_colors=gipr_tm_index_groups,
    )
    plt.title(f'GIPR co-mutation ratio when GLP1R is mutated')
    # Show the plot
    plt.savefig(script_dir / 'step5b_heatmap.png', dpi=300)
    plt.close()

    # only tm heatmap
    sns.set(font_scale=0.5)
    sns.set_style("whitegrid")
    # Create a heatmap
    sns.heatmap(
        df_only_tm,
        cmap="coolwarm",  # Color map for better visualization
        annot=False,  # Show cell values
        # fmt=".2f",  # Format cell values
        cbar=True,  # Show color bar
        square=True,  # Make cells square
        # deal with nans
        mask=np.isnan(df_only_tm),
        # color by group
        row_colors=glp1r_tm_index_groups,
        col_colors=gipr_tm_index_groups,
    )
    plt.title(f'GIPR co-mutation ratio when GLP1R is mutated (only TM residues)')
    # Show the plot
    plt.savefig(script_dir / 'step5b_heatmap_tm.png', dpi=300)


    # # create seabron clustermap
    # sns.set(font_scale=0.5)
    # sns.set_style("whitegrid")
    # # Create a clustermap
    # sns.clustermap(
    #     df,
    #     method="average",  # Clustering method (adjust if needed)
    #     metric="euclidean",  # Distance metric
    #     cmap="coolwarm",  # Color map for better visualization
    #     linewidths=0.5,  # Line width between cells
    #     figsize=(10, 10),  # Adjust figure size
    #     # standard_scale=1,  # Normalize each column
    #     # deal with nans
    #     mask=np.isnan(df),
    # )
    # plt.title(f'euclidean avg clustermap of GLP1R and GIPR mutations')
    # # Show the plot
    # plt.savefig(script_dir / 'step5b_clustermap.png', dpi=300)
    # plt.close()


    # # only tm clustermap
    # sns.set(font_scale=0.5)
    # sns.set_style("whitegrid")
    # # Create a clustermap
    # sns.clustermap(
    #     df_only_tm,
    #     method="average",  # Clustering method (adjust if needed)
    #     metric="euclidean",  # Distance metric
    #     cmap="coolwarm",  # Color map for better visualization
    #     linewidths=0.5,  # Line width between cells
    #     figsize=(10, 10),  # Adjust figure size
    #     standard_scale=1, # Normalize each column
    #     # deal with nans
    #     mask=np.isnan(df_only_tm),
    # )
    # plt.title(f'euclidean avg clustermap of GLP1R and GIPR mutations (only TM residues)')
    # # Show the plot
    # plt.savefig(script_dir / 'step5b_clustermap_tm.png', dpi=300)
    # plt.close()
