"""We have a lot of variant data for GLP1R and GIPR,
plot some overview of this data.

for example:
var_id	ref_aa3	pos_aa	alt_aa3	uniprot_id	gene	protein_variant	am_class	chr	pos	ref	alt	ref_aa1	alt_aa1	Counts	AC	AF	N_HOM	am_pathogenicity	PAI3D
Ala126Gly	Ala	126	Gly	P48546	GIPR	A126G	benign	chr19	45672947	C	G	A	G	4	6	6.11583389394e-6	0	0.0894	0.4005443243682385
Ala150Thr	Ala	150	Thr	P48546	GIPR	A150T	benign	chr19	45674137	G	A	A	T	222	2619	0.0026697792003913803	0	0.0829	0.32260909367352736
Ala185Gly	Ala	185	Gly	P48546	GIPR	A185G	benign	chr19	45674747	C	G	A	G	10	6	6.11655711228e-6	0	0.1286	0.4534999146312475
Ala185Val	Ala	185	Val	P48546	GIPR	A185V	benign	chr19	45674747	C	T	A	V	12	80	8.155409483109999e-5	0	0.0805	0.2649358920753002
Ala207Gly	Ala	207	Gly	P48546	GIPR	A207G	benign	chr19	45674813	C	G	A	G	14	128	1.3049107461432e-4	0	0.0798	0.3866692834347486


"""

import pathlib
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from settings import get_tm1_dict, get_tm4_dict, get_protein_lens


def load_variant_tsv(variant_tsv_p: pathlib.Path):
    """Load the variant tsv file into a pandas dataframe.

    variant.tsv has the following columns:
    'var_id', str, the ID. For example 'Ala126Gly'
    'ref_aa3', str, the reference amino acid. For example 'Ala'
    'pos_aa', int, the position of the amino acid. For example 126
    'alt_aa3', str, the alternative amino acid. For example 'Gly'
    'uniprot_id', str, the uniprot ID. For example 'P48546'
    'gene', str, the gene. For example 'GIPR' or 'GLP1R'
    'protein_variant', str, the protein variant. For example 'A126G'
    'am_class', str, the class of the variant predicted by AlphaMissense. For example 'benign'
    'chr', str, the chromosome. For example 'chr19'
    'pos', int, the position. For example 45672947
    'ref', str, the reference nucleotide. For example 'C'
    'alt', str, the alternative nucleotide. For example 'G'
    'ref_aa1', str, the reference amino acid. For example 'A'
    'alt_aa1', str, the alternative amino acid. For example 'G'
    'Counts', int, ?
    'AC', int, ?
    'AF', float, the allele frequency. For example 6.11583389394e-6
    'N_HOM', int, the number of homozygous individuals. For example 0
    'am_pathogenicity', float, the pathogenicity score predicted by AlphaMissense. For example 0.0894
    'PAI3D', float, the pathogenicity score predicted by Primate AI. For example 0.4005443243682385
    """
    df = pd.read_csv(
        variant_tsv_p,
        sep="\t",
        header=0,
    )
    # print head with types
    print(df.head())
    print(df.dtypes)

    # drop if nan at pos_aa
    dropped_rows = df[df["pos_aa"].isna()]
    print(f"Dropped {len(dropped_rows)} rows with nan at pos_aa")
    print(dropped_rows)

    df = df.dropna(subset=["pos_aa"])
    df["pos_aa"] = df["pos_aa"].astype(int, errors="ignore")
    df["pos"] = df["pos"].astype(int, errors="ignore")
    df["Counts"] = df["Counts"].astype(int, errors="ignore")
    df["AC"] = df["AC"].astype(int, errors="ignore")
    df["N_HOM"] = df["N_HOM"].astype(int, errors="ignore")
    return df


def plot_variant_counts(df: pd.DataFrame, plot_path: pathlib.Path):
    """Plot the counts of the variants for each gene."""
    # group the data by gene
    gene_counts = df.groupby("gene")["var_id"].count()
    gene_counts = gene_counts.sort_values(ascending=False)
    gene_counts.plot(kind="bar")
    plt.title("Number of variants per gene")
    plt.xlabel("Gene")
    plt.ylabel("Number of variants")
    plt.tight_layout()
    plt.savefig(plot_path)
    plt.close()


def plot_variant_positions_and_frequency(
    plot_df: pd.DataFrame,
    plot_path: pathlib.Path,
    color_by_column: str = "am_pathogenicity",
):
    """Plot the frequency of the variants and where they are in the protein."""
    # copy df
    df = plot_df.copy()

    # get unique genes
    genes = df["gene"].unique()
    # protein lengths
    protein_lens = get_protein_lens()
    # AF2 observation: the heterodimer with glp1 forms through TM4-TM1 pairing. While with the gip peptide, this interface does not exist and TM4 is turned away
    tm1_dict = get_tm1_dict()
    tm4_dict = get_tm4_dict()

    # create a figure with subplots, one for each gene
    fig, axes = plt.subplots(1, len(genes), figsize=(15, 7))

    # get min max vals for color bar
    min_color_by = df[color_by_column].min()
    max_color_by = df[color_by_column].max()

    # change AF to log scale
    df["AF"] = np.log10(df["AF"])
    for gene, ax in zip(genes, axes):
        # get the variants for this gene
        gene_df = df[df["gene"] == gene]
        # get the allele frequencies
        af = gene_df["AF"]
        # get the positions
        pos = gene_df["pos_aa"]
        # get the colors
        vals_to_norm = gene_df[color_by_column]
        cmap = plt.colormaps[
            "plasma"
        ]  # other cmaps: "plasma", "inferno", "magma", "cividis"
        colors = cmap((vals_to_norm - min_color_by) / (max_color_by - min_color_by))

        # highlight tm1 and tm4, depending on the gene
        if gene == "GIPR":
            tm1 = tm1_dict[gene]
            ax.axvspan(tm1[0], tm1[1], color="lightgrey", alpha=0.5)
            ax.text(tm1[0] + 3, max(af) - 0.2, "TM1", fontsize=8)
        elif gene == "GLP1R":
            tm4 = tm4_dict[gene]
            ax.axvspan(tm4[0], tm4[1], color="lightgrey", alpha=0.5)
            ax.text(tm4[0] + 3, max(af) - 0.2, "TM4", fontsize=8)
        else:
            raise ValueError(f"Unknown gene {gene}")

        # plot the allele frequencies
        # summarize af
        ax.scatter(pos, af, s=5, c=colors)
        ax.set_title(gene)
        ax.set_xlabel("Position")
        ax.set_ylabel("log10 Allele Frequency")
        ax.set_xlim(0, protein_lens[gene])

        # give xticks every 50
        ax.set_xticks(np.arange(0, protein_lens[gene] + 1, 50))

    # add a color bar
    sm = plt.cm.ScalarMappable(
        cmap=cmap, norm=plt.Normalize(vmin=min_color_by, vmax=max_color_by)
    )
    cbar = plt.colorbar(
        sm,
        ax=axes.ravel().tolist(),
        orientation="horizontal",
        # place below the subplots
        pad=0.2,
        # scale the color bar
        shrink=0.3,
    )
    # add a title
    plt.suptitle(
        f"GLP1R and GIPR variant (UK BioBank) frequency and predicted {color_by_column}"
    )
    cbar.set_label(color_by_column)
    cbar.set_ticks([min_color_by, max_color_by])
    plt.savefig(plot_path, dpi=300)
    plt.close()


def plot_correlation_PAI3D_am_pathogenicity(df: pd.DataFrame, plot_path: pathlib.Path):
    """Plot the correlation between the PAI3D and am_pathogenicity scores."""
    sns.jointplot(x="PAI3D", y="am_pathogenicity", data=df, kind="reg")
    plt.savefig(plot_path)
    plt.close()


if __name__ == "__main__":
    VARIANT_TSV_P = pathlib.Path("data/variants.tsv")
    VARIANT_DF = load_variant_tsv(VARIANT_TSV_P)

    PLOT_DIR = pathlib.Path("plots/variants/")
    if not PLOT_DIR.exists():
        PLOT_DIR.mkdir()

    # set matplotlib style for publication
    print(plt.style.available)
    plt.style.use("seaborn-v0_8-ticks")

    plot_variant_counts(VARIANT_DF, PLOT_DIR / "variant_counts.png")

    plot_variant_positions_and_frequency(
        VARIANT_DF,
        PLOT_DIR / "variant_positions_and_frequency_PAI3D.png",
        color_by_column="PAI3D",
    )

    plot_variant_positions_and_frequency(
        VARIANT_DF,
        PLOT_DIR / "variant_positions_and_frequency_am_pathogenicity.png",
        color_by_column="am_pathogenicity",
    )

    plot_correlation_PAI3D_am_pathogenicity(
        VARIANT_DF, PLOT_DIR / "correlation_PAI3D_am_pathogenicity.png"
    )
    # separate plots for each gene
    for GENE in VARIANT_DF["gene"].unique():
        GENE_DF = VARIANT_DF[VARIANT_DF["gene"] == GENE]
        plot_correlation_PAI3D_am_pathogenicity(
            GENE_DF, PLOT_DIR / f"correlation_PAI3D_am_pathogenicity_{GENE}.png"
        )
