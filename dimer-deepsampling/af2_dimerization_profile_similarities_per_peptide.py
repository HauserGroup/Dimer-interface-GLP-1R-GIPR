"""
This script considers the research question:

How similar are the AF2 dimerization profiles between the different peptides?
For example, is GLP-1 dimerization more similar to Exendin-4 than GIP?


Process;
    - the dimer distances (angstrom between CA) have already been calculated
    - for each peptide, seed, model combination, take the best (highest iptm) sample
    - calculate the euclidean distance between the dimer distance maps of each peptide pair
    - plot boxplots for each peptide, and the 50 distance values per other peptide

"""
import os
import pathlib
import pandas as pd
from tmtools import tm_align

import numpy as np
from tqdm import tqdm
from scipy.spatial.distance import pdist, squareform

import seaborn as sns
import matplotlib.pyplot as plt

import hdbscan
import umap
from sklearn.preprocessing import StandardScaler

from scipy.spatial.distance import pdist, squareform
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from statannotations.Annotator import Annotator
import matplotlib as mpl

# np set seed
np.random.seed(1)

# svg editing
plt.rcParams['svg.fonttype'] = 'none'
# fontsize
plt.rcParams.update({'font.size': 8})
# style
sns.set_style("whitegrid")

COLORS = {
    "Exendin-4" : "#5cbb60",
    "Exendin-9" : "#ffffff",
    "GIP" : "#ed903e",
    "GLP-1" : "#79a1d4",
    "Semaglutide" : "#3548a1",
    "Tirzepatide" : "#e53d35",
}



def plot_boxplot_dimerization_profile_distances(plot_df: pd.DataFrame,
                                                out_p: pathlib.Path,
                                                palette : dict):
    plt.figure(figsize=(10, 4))
    ax = sns.boxplot(
        data=plot_df, 
        x='peptide', 
        y='euclidean_distance', 
        hue='other_peptide',
        palette=palette,
    )
    # add stripplot
    sns.stripplot(
        data=plot_df,
        x='peptide',
        y='euclidean_distance',
        hue='other_peptide',
        palette={x : "black" for x in palette.keys()},
        dodge=True,
        size=1,
        alpha=0.3,
        ax=ax
    )

    ax.set_title("Peptide Dimerization Profile Euclidean Distances")
    ax.set_ylabel("Euclidean Distance")
    ax.set_xlabel("Peptide")
    plt.legend(title='Other Peptide', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(out_p, dpi=300)
    plt.close()


# alpha-carbon distance -----------------------------------------------------------------
def get_euclidean_distance_matrix(cmap1: pd.DataFrame, cmap2: pd.DataFrame) -> float:
    # assert the keys and cols match
    assert (cmap1.index == cmap2.index).all(), "CMap row indices do not match"
    assert (cmap1.columns == cmap2.columns).all(), "CMap column indices do not match"

    # calculate euclidean distance
    diff = cmap1.values - cmap2.values
    euclidean_distance = np.sqrt(np.sum(diff ** 2))
    return euclidean_distance


def calculate_peptide_dimerization_profile_distances(df: pd.DataFrame, out_p: pathlib.Path):
    # sanity
    assert df['id'].nunique() == len(df), "IDs are not unique per row"

    # keep a stable order
    df2 = df[['id', 'dimer_distances']].drop_duplicates('id').sort_values('id')
    ids = df2['id'].to_numpy()

    # 1) Load each map ONCE, convert to numpy, flatten
    vectors = []
    shape0 = None
    for path in tqdm(df2['dimer_distances'], desc="Loading distance maps"):
        arr = pd.read_csv(path, index_col=0).to_numpy(dtype=np.float32, copy=False)
        if shape0 is None:
            shape0 = arr.shape
        elif arr.shape != shape0:
            raise ValueError(f"Shape mismatch: expected {shape0}, got {arr.shape} for {path}")

        vectors.append(arr.ravel())

    # X: (n_samples, n_features)
    X = np.vstack(vectors)  # float32 to reduce memory + often faster

    # 2) Compute all pairwise Euclidean distances in optimized C
    dists = squareform(pdist(X, metric="euclidean")).astype(np.float32, copy=False)
    np.fill_diagonal(dists, 0.0)

    # 3) Save
    distances_df = pd.DataFrame(dists, index=ids, columns=ids)
    distances_df.to_csv(out_p)
# alpha-carbon distance -----------------------------------------------------------------



# rmsd ----------------------------------------------------------------------------------    
from tmtools import tm_align
from tmtools.io import get_structure, get_residue_data  # needs biopython

def _load_coords_seq(pdb_p: str | pathlib.Path, chain_id: str | None = None):
    """Load (coords, seq) once from a PDB file (optionally choose a chain)."""
    pdb_p = pathlib.Path(pdb_p)
    s = get_structure(str(pdb_p))

    if chain_id is None:
        chain = next(s.get_chains())
    else:
        chain = None
        for ch in s.get_chains():
            if ch.id == chain_id:
                chain = ch
                break
        if chain is None:
            raise ValueError(f"Chain {chain_id!r} not found in {pdb_p}")

    coords, seq = get_residue_data(chain)
    return coords, seq


def _calc_rmsd_from_coords(coords1: np.ndarray, seq1: str, coords2: np.ndarray, seq2: str) -> float:
    res = tm_align(coords1, coords2, seq1, seq2)
    return float(res.rmsd)


def calculate_peptide_dimerization_profile_rmsd_distances(
    df: pd.DataFrame,
    out_p: pathlib.Path,
    chain_id: str | None = None,  # set if you want e.g. "A"
):
    # sanity
    assert df["id"].nunique() == len(df), "IDs are not unique per row"

    # keep a stable order
    df2 = df[["id", "pdb_p"]].drop_duplicates("id").sort_values("id").reset_index(drop=True)
    ids = df2["id"].to_numpy()
    n = len(ids)

    # ---- parse once (this is where apply actually helps) ----
    # This avoids re-reading PDBs N*(N-1)/2 times.
    data = df2["pdb_p"].apply(lambda p: _load_coords_seq(p, chain_id=chain_id)).tolist()
    # data[i] = (coords_i, seq_i)

    rmsd_matrix = np.zeros((n, n), dtype=np.float32)

    # ---- pairwise TM-align on cached coords ----
    for i in tqdm(range(n), desc="Calculating pairwise RMSD"):
        coords1, seq1 = data[i]
        for j in range(i + 1, n):
            coords2, seq2 = data[j]
            rmsd = _calc_rmsd_from_coords(coords1, seq1, coords2, seq2)
            rmsd_matrix[i, j] = rmsd
            rmsd_matrix[j, i] = rmsd  # symmetric

    rmsd_df = pd.DataFrame(rmsd_matrix, index=ids, columns=ids)
    rmsd_df.to_csv(out_p)
# rmsd ----------------------------------------------------------------------------------


def _fill_in_plot_df(euclid_distance_df, 
                     df: pd.DataFrame, 
                     peptide: str, 
                     other_peptide: str):
    """
    df decides which ids are included
    """
    plot_df = []
    ids = df.loc[df['peptide'] == peptide, 'id']
    other_ids = df.loc[df['peptide'] == other_peptide, 'id']
    # get distances
    for id1 in ids:
        for id2 in other_ids:
            dist = euclid_distance_df.at[id1, id2]
            plot_df.append({
                'id1' : id1,
                'id2' : id2,
                'peptide': peptide,
                'other_peptide': other_peptide,
                'euclidean_distance': dist
            })
    return plot_df



# plot --------------------------------------------------------------------------------
def plot_sorted_median_boxplots_per_peptide(
    plot_df: pd.DataFrame,
    out_p,
    palette: dict,
    ncols: int = 3,
    show_points: bool = True,
    point_size: float = 1.2,
    point_alpha: float = 0.25,
    box_width: float = 0.65,
):
    """
    One subplot per peptide. Within each subplot, boxplots for `other_peptide`
    are sorted by the median euclidean distance (ascending).

    Expects columns:
        - peptide
        - other_peptide
        - euclidean_distance
    """

    # --- Nature-ish style tweaks (local, not global) ---
    sns.set_theme(style="white", context="paper")
    plt.rcParams.update({
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "legend.frameon": False,
        "svg.fonttype": "none",  # keep editable text in svg
    })

    # Ensure deterministic peptide order (optional: keep original appearance order)
    peptides = list(pd.unique(plot_df["peptide"]))

    n = len(peptides)
    ncols = max(1, min(ncols, n))
    nrows = int(np.ceil(n / ncols))

    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(3.2 * ncols, 2.6 * nrows),
        sharey=True,
        constrained_layout=True
    )
    axes = np.atleast_1d(axes).ravel()

    # Global y-limits with a small margin for consistency
    y = plot_df["euclidean_distance"].to_numpy()
    y_min, y_max = np.nanmin(y), np.nanmax(y)
    pad = 0.03 * (y_max - y_min if y_max > y_min else 1.0)
    ylims = (y_min - pad, y_max + pad)

    for i, pep in enumerate(peptides):
        ax = axes[i]
        sub = plot_df.loc[plot_df["peptide"] == pep].copy()

        # Compute medians per other_peptide and sort
        med = (
            sub.groupby("other_peptide", as_index=False)["euclidean_distance"]
            .median()
            .sort_values("euclidean_distance", ascending=True)
        )
        order = med["other_peptide"].tolist()

        # Build a palette only for the categories present in this subplot
        pal = {k: palette.get(k, "#666666") for k in order}

        # Boxplot: thin lines, no fliers, median line visible
        sns.boxplot(
            data=sub,
            x="other_peptide",
            y="euclidean_distance",
            order=order,
            palette=pal,
            hue="other_peptide",
            width=box_width,
            fliersize=0,
            linewidth=0.9,
            ax=ax,
            medianprops=dict(color="black", linewidth=1.2),
            whiskerprops=dict(linewidth=0.9),
            capprops=dict(linewidth=0.9),
            boxprops=dict(linewidth=0.9),
        )

        # Optional points (small + subtle), jittered
        if show_points:
            sns.stripplot(
                data=sub,
                x="other_peptide",
                y="euclidean_distance",
                order=order,
                color="black",
                size=point_size,
                alpha=point_alpha,
                jitter=0.18,
                ax=ax
            )

        # Title + axes formatting
        ax.set_title(pep, pad=3,
                     # add border
                     bbox=dict(facecolor=palette.get(pep, "#ffffff"), edgecolor='black', boxstyle='round,pad=0.2')
        )
                     
        ax.set_xlabel("")  # cleaner: categories explain themselves
        if i % ncols == 0:
            ax.set_ylabel("Euclidean distance (Å)")
        else:
            ax.set_ylabel("")

        ax.set_ylim(*ylims)

        # Clean spines (Nature-ish)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        # Subtle horizontal grid only
        ax.grid(True, axis="y", linewidth=0.5, alpha=0.25)
        ax.grid(False, axis="x")

        # Rotate category labels for readability
        ax.tick_params(axis="x", rotation=45)
        for tick in ax.get_xticklabels():
            tick.set_ha("right")

    # Turn off any unused axes
    for j in range(n, len(axes)):
        axes[j].axis("off")

    fig.suptitle("Peptide dimerization profile distances", y=1.02, fontsize=10)
    fig.savefig(out_p, dpi=600, bbox_inches="tight")
    plt.close(fig)

MM_PER_INCH = 25.4
def mm_to_in(mm): 
    return mm / MM_PER_INCH


def plot_dimerization_cumalitive_euclidean_distance(
    plot_df: pd.DataFrame,
    out_p,
    palette: dict,
    ncols: int = 3,
    show_points: bool = False,
):
    """
    For each peptide (subplot), plot empirical CDF curves of euclidean_distance
    for each other_peptide:
        y = fraction of distances <= x

    Expects columns:
        - peptide
        - other_peptide
        - euclidean_distance
    """
    required = {"peptide", "other_peptide", "euclidean_distance"}
    missing = required - set(plot_df.columns)
    if missing:
        raise ValueError(f"plot_df missing columns: {missing}")

    # --- Supplementary-figure style (clean, print-friendly) ---
    sns.set_theme(style="white", context="paper")
    plt.rcParams.update({
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,

        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.major.size": 3,
        "ytick.major.size": 3,

        "legend.frameon": False,
        "svg.fonttype": "none",
        # makes joins/caps look cleaner in vector and raster export
        "lines.solid_joinstyle": "round",
        "lines.solid_capstyle": "round",
    })

    peptides = list(pd.unique(plot_df["peptide"]))
    n = len(peptides)
    ncols = max(1, min(ncols, n))
    nrows = int(np.ceil(n / ncols))

    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(3.2 * ncols, 2.6 * nrows),
        sharex=True,
        sharey=True,
        constrained_layout=True
    )
    axes = np.atleast_1d(axes).ravel()

    # Global x limits (consistent across panels)
    x = plot_df["euclidean_distance"].to_numpy()
    x = x[np.isfinite(x)]
    x_min, x_max = np.min(x), np.max(x)
    pad = 0.03 * (x_max - x_min if x_max > x_min else 1.0)
    xlims = (x_min - pad, x_max + pad)

    # Helper: ECDF points as a true step function (publication standard for CDFs)
    def _ecdf_step(arr: np.ndarray):
        arr = np.asarray(arr, dtype=float)
        arr = arr[np.isfinite(arr)]
        if arr.size == 0:
            return np.array([]), np.array([])
        xs = np.sort(arr)
        ys = np.arange(1, xs.size + 1, dtype=float) / xs.size
        return xs, ys

    for i, pep in enumerate(peptides):
        ax = axes[i]
        sub = plot_df.loc[plot_df["peptide"] == pep].copy()

        # Sort curves by median distance (low -> high) so ordering is meaningful
        med = (
            sub.groupby("other_peptide", as_index=False)["euclidean_distance"]
            .median()
            .sort_values("euclidean_distance", ascending=True)
        )
        order = med["other_peptide"].tolist()

        for other in order:
            arr = sub.loc[sub["other_peptide"] == other, "euclidean_distance"].to_numpy()
            xs, ys = _ecdf_step(arr)
            if xs.size == 0:
                continue

            color = palette.get(other, "#666666")

            # High-quality ECDF rendering:
            # - use a step plot (empirical CDF)
            # - slightly thicker but not heavy
            # - round caps/joins (rcParams) for cleaner appearance
            ax.step(
                xs, ys,
                where="post",
                linewidth=1.25,
                color=color,
                alpha=0.95,
                label=other,
                zorder=2,
                # style tweaks
                solid_capstyle="round",
                solid_joinstyle="round",

            )

            if show_points:
                # optional: *very* subtle points, mainly for debugging
                ax.plot(
                    xs, ys,
                    linestyle="None",
                    marker="o",
                    markersize=1.0,
                    alpha=0.18,
                    color=color,
                    zorder=1
                )

        ax.set_title(
            pep,
            pad=6,
            bbox=dict(
                facecolor=palette.get(pep, "#ffffff"),
                edgecolor="black",
                linewidth=0.6,
                boxstyle="round,pad=0.22",
            ),
        )
        ax.set_xlim(*xlims)
        ax.set_ylim(0, 1.0)

        # Clean spines (keep left/bottom for print)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        # Ticks: clean and consistent
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.tick_params(axis="both", which="major", direction="out")

        if i % ncols == 0:
            ax.set_ylabel("Fraction ≤ distance")
        else:
            ax.set_ylabel("")

        if i >= (nrows - 1) * ncols:
            ax.set_xlabel("Euclidean distance (Å)")
        else:
            ax.set_xlabel("")

        # Compact legend inside panel; keep line samples long enough to see thickness
        ax.legend(
            title=None,
            fontsize=6.5,
            loc="lower right",
            frameon=False,
            handlelength=2.4,
            handletextpad=0.6,
            borderaxespad=0.25
        )

        ax.grid(True, axis="both", linewidth=1.0, alpha=0.25,
                color="gray", zorder=0,
                linestyle="--")

    # Turn off unused axes
    for j in range(n, len(axes)):
        axes[j].axis("off")

    fig.suptitle("Empirical CDF of peptide dimerization profile distances", y=1.02, fontsize=10)
    fig.savefig(out_p, dpi=600, bbox_inches="tight")
    plt.close(fig)
# plot --------------------------------------------------------------------------------



# residue participation analysis ------------------------------------------------------
def plot_residue_participation_analysis_per_peptide(
    df: pd.DataFrame,
    out_dir: pathlib.Path,
    palette: dict,
    angstrom_cutoff: float = 8.0,
    ncols: int = 3,
    ):
    
    # instead of angstrom distances, use angstrom cutoff and label contacts 1/0
    def _get_participation_array(dimer_distances_p):
        cmap = pd.read_csv(dimer_distances_p, index_col=0)
        # binarize based on cutoff
        contact_map = (cmap.values <= angstrom_cutoff).astype(int)
        # sum over rows and cols to get participation per residue.
        # 463x466 -> 463+466 = 929 array
        participation_array = np.concatenate([
            contact_map.sum(axis=0),  # chain A residues
            contact_map.sum(axis=1),  # chain B residues
        ])
        return participation_array.tolist()
    
    # save df
    df['participation_array'] = df['dimer_distances'].apply(_get_participation_array)
    df.to_csv(out_dir / "residue_participation_arrays_per_model.csv", index=False)


    # run HDBScan and color on peptide
    participation_matrix = np.vstack(df['participation_array'].to_numpy())
    scaler = StandardScaler()
    embedder = umap.UMAP(random_state=2026,
                         n_neighbors=25,
                         min_dist=0.25,
                         metric='euclidean')
    participation_matrix_scaled = scaler.fit_transform(participation_matrix)
    participation_embedding = embedder.fit_transform(participation_matrix_scaled)
    df_embedding = pd.DataFrame(participation_embedding, columns=['UMAP1', 'UMAP2'])
    df_embedding = pd.concat([df.reset_index(drop=True), df_embedding], axis=1)

    # TODO: use same umap for ref
    ref_p = pathlib.Path(__file__).parent.resolve() / "reference/reference_distance_map_A_B.csv"
    ref_participation_array = _get_participation_array(ref_p)
    print("Reference participation array loaded from", ref_p)
    ref_participation_array_scaled = scaler.transform([ref_participation_array])
    ref_embedding = embedder.transform(ref_participation_array_scaled)
    print("Reference participation array UMAP embedding:", ref_embedding)
    
    # plot
    plt.figure(figsize=(6, 5))
    ax = sns.scatterplot(
        data=df_embedding,
        x='UMAP1',
        y='UMAP2',
        hue='peptide',
        palette=palette,
        alpha=0.7,
        s=20,
        edgecolor=None,
    )
    ax.set_title("Residue Participation Analysis of Peptide Dimerization Profiles")
    plt.legend(title='Peptide', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(out_dir / "residue_participation_analysis_per_peptide.png", dpi=300)
    plt.close()
 
    # theme/rc MUST be set before plotting
    sns.set_theme(style="white", context="paper")
    plt.rcParams.update({
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 8,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "legend.frameon": False,
        "svg.fonttype": "none",
        "lines.solid_joinstyle": "round",
        "lines.solid_capstyle": "round",
    })


    # ---------------------------------------------------------------------
    # -------- global limits (start at 0, add small headroom) --------
    x = df_embedding["UMAP1"].to_numpy()
    y = df_embedding["UMAP2"].to_numpy()

    xmin = float(np.nanmin(x))
    ymin = float(np.nanmin(y))
    xmax = float(np.nanmax(x))
    ymax = float(np.nanmax(y))

    # padding as fraction of range
    pad = 0.5
    xr = max(xmax - xmin, 1e-9)
    yr = max(ymax - ymin, 1e-9)
    xlim = (xmin - pad * xr, xmax + pad * xr)
    ylim = (ymin - pad * yr, ymax + pad * yr)


    # -------- FacetGrid --------
    facet_h_mm = 210
    facet_w_mm = 50  # (note: 50mm is very narrow; equal-aspect UMAPs will look squeezed)

    g = sns.FacetGrid(
        df_embedding,
        col="peptide",
        col_wrap=ncols,
        sharex=True,
        sharey=True,
        height=1.5,
        aspect=1.0,
    )

    # IMPORTANT: set limits BEFORE kdeplot so each facet is drawn in the same coord system
    for ax in g.axes.flatten():
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)

    # Draw KDE with a stable grid and no extrapolation outside data range
    # - cut=0 prevents KDE extending beyond data range
    # - gridsize stabilizes the evaluation mesh so plots look consistent facet-to-facet
    g.map_dataframe(
        sns.kdeplot,
        x="UMAP1", 
        y="UMAP2",
        fill=True,
        levels=10,
        thresh=0.10,
        cmap="mako",
        #cut=0,
        #gridsize=200,
    )

    # Reference point (will now align with the same locked axes)
    for ax in g.axes.flatten():
        ax.plot(
            ref_embedding[0, 0], ref_embedding[0, 1],
            marker="X",
            color="white",
            markersize=6,
            markeredgecolor="black",
            markeredgewidth=0.7,
            zorder=5,
        )

    # Titles
    for ax, title in zip(g.axes.flatten(), g.col_names):
        ax.set_title(
            title,
            pad=2,
            bbox=dict(
                facecolor=palette.get(title, "#ffffff"),
                edgecolor="black",
                linewidth=0.6,
                boxstyle="round,pad=0.18",
            ),
        )

    # Keep x ticks/labels on EVERY facet (seaborn hides them by default with sharex)
    for ax in g.axes.flatten():
        ax.tick_params(axis="x", which="both", bottom=True, labelbottom=True, pad=1)
        for lab in ax.get_xticklabels():
            lab.set_visible(True)

    # Layout
    g.fig.subplots_adjust(left=0.06, right=0.90, bottom=0.08, top=0.90, wspace=0.15, hspace=0.18)

    # Shared colorbar (relative density label only; KDE itself isn't normalized to [0,1] here)
    norm = mpl.colors.Normalize(vmin=0, vmax=1)
    sm = mpl.cm.ScalarMappable(cmap="mako", norm=norm)
    sm.set_array([])
    cbar = g.fig.colorbar(sm, ax=g.axes, fraction=0.018, pad=0.01)
    cbar.set_label("Relative density (within peptide)", fontsize=7, labelpad=2)
    cbar.ax.tick_params(labelsize=6, pad=1, length=2)

    # Suptitle
    g.fig.suptitle("Density of Residue Participation Profiles per Peptide", y=0.97, fontsize=9)

    # add grid with alpha=0.25
    for ax in g.axes.flatten():
        ax.grid(True, axis="both", linewidth=0.5, alpha=0.25,
                color="gray", zorder=0,
                linestyle="--")

    # Save
    g.fig.savefig(out_dir / "residue_participation_umap_density_faceted.png", dpi=300, bbox_inches="tight")
    g.fig.savefig(out_dir / "residue_participation_umap_density_faceted.svg", dpi=300, bbox_inches="tight")
    plt.close(g.fig)

    # ---------------------------------------------------------------------

    # # ###
    # # _plot_participation_ecf(
    # #     df_embedding=df_embedding,
    # #     ref_participation_array=ref_participation_array,
    # #     out_p=out_dir / "residue_participation_ecdf_per_peptide.png",
    # #     palette=palette,
    # #     ncols=3,
    # #     show_points=False,
    # #     out_dir=out_dir,
    # # )
    

    #call jaccard
    jaccard_sim = _run_participation_jaccard_similarity(
        participation_matrix=participation_matrix,
        ids=df["id"].to_numpy(),   # <-- critical
        out_p=out_dir / "residue_participation_jaccard_similarity.csv",
    )
    longform_df = build_per_sample_mean_jaccard(
        jaccard_sim=jaccard_sim,
        df_meta=df,
        value_col="jaccard_mean",
    )
    longform_df.to_csv(
        out_dir / "residue_participation_jaccard_similarity_longform.csv",
        index=False,
    )
    plot_within_vs_other_per_peptide_aggregated(
        per_sample_df=longform_df,
        out_p=out_dir / "residue_participation_jaccard_within_vs_other_boxplots.png",
        palette=palette,
        ncols=1,
        show_points=True,
        point_size=1.5,
        point_alpha=0.25,
        stat_test="Mann-Whitney",
        # bonferroni correction for multiple comparisons
        p_adjust="bonferroni",
        text_format="star",
        loc="inside",
        value_col="jaccard_mean",
        min_n=3,
    )

def _run_participation_jaccard_similarity(
    participation_matrix: np.ndarray,
    ids: np.ndarray | list,
    out_p: pathlib.Path,
):
    """
    Compute all-vs-all Jaccard similarity between binarized participation profiles.
    Uses provided `ids` for matrix index/columns (must match row order of participation_matrix).
    """
    out_p = pathlib.Path(out_p)
    out_p.parent.mkdir(parents=True, exist_ok=True)

    X = np.asarray(participation_matrix)
    if X.ndim != 2:
        raise ValueError(f"participation_matrix must be 2D, got shape {X.shape}")
    if len(ids) != X.shape[0]:
        raise ValueError(f"len(ids) ({len(ids)}) must equal n_samples ({X.shape[0]})")

    # Binarize counts -> presence/absence
    X_bin = (X > 0).astype(bool, copy=False)

    # pdist('jaccard') returns distance = 1 - similarity
    d = pdist(X_bin, metric="jaccard")
    D = squareform(d)
    np.fill_diagonal(D, 0.0)

    S = 1.0 - D
    np.fill_diagonal(S, 1.0)

    sim_df = pd.DataFrame(S, index=ids, columns=ids)
    sim_df.to_csv(out_p)
    return sim_df

from scipy.stats import wilcoxon

def build_per_sample_mean_jaccard(
    jaccard_sim: pd.DataFrame,
    df_meta: pd.DataFrame,
    value_col: str = "jaccard_mean",
) -> pd.DataFrame:
    """
    Build an aggregated dataframe where each AF2 structure (sample) contributes ONE value
    per comparison group for its focal peptide, using the MEAN (average) similarity.

    Output columns:
      focal         : focal peptide A
      other         : 'within' or peptide B
      id            : sample id (from df_meta['id'])
      jaccard_mean  : mean similarity of that sample to group 'other'
    """
    # Align metadata to similarity matrix ordering
    ids = jaccard_sim.index.to_numpy()
    meta = df_meta.set_index("id").loc[ids].reset_index()  # preserves order
    peptides = meta["peptide"].to_numpy()

    S = jaccard_sim.to_numpy(dtype=float, copy=False)
    n = S.shape[0]
    if S.shape[1] != n:
        raise ValueError("jaccard_sim must be square")

    unique_peps = sorted(pd.unique(peptides))
    idx_by_pep = {p: np.flatnonzero(peptides == p) for p in unique_peps}

    rows = []
    for A in unique_peps:
        idxA = idx_by_pep[A]
        if idxA.size == 0:
            continue

        # --- within(A): for each i in A, mean to A excluding self ---
        for i in idxA:
            othersA = idxA[idxA != i]
            if othersA.size == 0:
                continue
            val = np.nanmean(S[i, othersA])
            rows.append((A, "within", ids[i], float(val)))

        # --- A vs B: for each i in A, mean to all samples in B ---
        for B in unique_peps:
            if B == A:
                continue
            idxB = idx_by_pep[B]
            if idxB.size == 0:
                continue
            for i in idxA:
                val = np.nanmean(S[i, idxB])
                rows.append((A, B, ids[i], float(val)))

    out = pd.DataFrame(rows, columns=["focal", "other", "id", value_col])
    return out



def plot_within_vs_other_per_peptide_aggregated(
    per_sample_df: pd.DataFrame,
    out_p: pathlib.Path,
    palette: dict,
    ncols: int = 1,
    show_points: bool = True,
    point_size: float = 1.0,
    point_alpha: float = 0.25,
    stat_test: str = "Mann-Whitney",   # or try "t-test_welch" / "t-test_ind"
    p_adjust: str = "BH",              # "bonferroni", "holm", "BH", etc.
    text_format: str = "star",
    loc: str = "inside",
    value_col: str = "jaccard_median",
    min_n: int = 3,                    # skip stats if either group has < min_n samples
):
    """
    Same visual goal as your 1-vs-5 boxplots, but using per-sample aggregated values
    (one value per structure per comparison group), so the stats are not pseudo-replicated.

    Stats added: within vs each other peptide (within, B) per subplot using statannotations.
    """
    plot_df = per_sample_df.copy()
    if value_col not in plot_df.columns:
        raise ValueError(f"per_sample_df missing {value_col!r}")
    if not {"focal", "other"}.issubset(plot_df.columns):
        raise ValueError("per_sample_df must contain columns: focal, other")

    peptides = sorted(pd.unique(plot_df["focal"]))

    # --- style (match your paper-ish look) ---
    sns.set_theme(style="white", context="paper")
    plt.rcParams.update({
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "legend.frameon": False,
        "svg.fonttype": "none",
    })

    n = len(peptides)
    ncols_eff = max(1, min(ncols, n))
    nrows = int(np.ceil(n / ncols_eff))

    FIG_SCALE = 0.7  # ~50% area reduction

    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols_eff,
        figsize=(3.4 * ncols_eff * FIG_SCALE,
                2.5 * nrows * FIG_SCALE),
        sharey=True,
        constrained_layout=True,
    )
    axes = np.atleast_1d(axes).ravel()

    # global y-lims
    y = plot_df[value_col].to_numpy()
    y = y[np.isfinite(y)]
    y_min, y_max = (0.0, 1.0) if y.size == 0 else (
        max(0.0, float(np.min(y) - 0.02)),
        min(1.0, float(np.max(y) + 0.02)),
    )

    for i, A in enumerate(peptides):
        ax = axes[i]
        sub = plot_df[plot_df["focal"] == A].copy()
        if sub.empty:
            ax.axis("off")
            continue

        # Order: within first, then others by median similarity (high -> low)
        med = sub.groupby("other", as_index=False)[value_col].median()
        med_other = med[med["other"] != "within"].sort_values(value_col, ascending=False)
        order = (["within"] if (sub["other"] == "within").any() else []) + med_other["other"].tolist()

        pal = {lab: (palette.get(A, "#666666") if lab == "within" else palette.get(lab, "#666666"))
               for lab in order}

        sns.boxplot(
            data=sub,
            x="other",
            y=value_col,
            order=order,
            palette=pal,
            width=0.65,
            fliersize=0,
            linewidth=0.9,
            ax=ax,
            medianprops=dict(color="black", linewidth=1.2),
            whiskerprops=dict(linewidth=0.9),
            capprops=dict(linewidth=0.9),
            boxprops=dict(linewidth=0.9),
        )

        if show_points:
            sns.stripplot(
                data=sub,
                x="other",
                y=value_col,
                order=order,
                color="black",
                size=point_size,
                alpha=point_alpha,
                jitter=0.18,
                ax=ax,
            )

        # ---- stats: within vs each other peptide (1 vs 5) ----
        if "within" in order and len(order) > 1:
            pairs = [("within", lab) for lab in order if lab != "within"]

            good_pairs = []
            for a, b in pairs:
                na = sub.loc[sub["other"] == a, value_col].notna().sum()
                nb = sub.loc[sub["other"] == b, value_col].notna().sum()
                if na >= min_n and nb >= min_n:
                    good_pairs.append((a, b))

            if good_pairs:
                annot = Annotator(
                    ax=ax,
                    pairs=good_pairs,
                    data=sub,
                    x="other",
                    y=value_col,
                    order=order,
                )
                annot.configure(
                    test=stat_test,
                    comparisons_correction=p_adjust,  # correction within this subplot
                    text_format=text_format,
                    loc=loc,
                    verbose=0,
                )
                annot.apply_and_annotate()

        ax.set_title(
            A,
            pad=4,
            bbox=dict(
                facecolor=palette.get(A, "#ffffff"),
                edgecolor="black",
                linewidth=0.6,
                boxstyle="round,pad=0.22",
            ),
        )
        ax.set_ylim(y_min, y_max)
        ax.set_xlabel("")
        ax.set_ylabel("Mean Jaccard similarity" if i % ncols_eff == 0 else "")

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(True, axis="y", linewidth=0.5, alpha=0.25)
        ax.grid(False, axis="x")

        ax.tick_params(axis="x", rotation=45)
        for tick in ax.get_xticklabels():
            tick.set_ha("right")

        # make xticklabels sit closer to the axis
        ax.tick_params(axis="x", rotation=45, pad=0.5)   # try 0, 0.5, 1
        for tick in ax.get_xticklabels():
            tick.set_ha("right")
            tick.set_va("top")  # helps when rotated; keeps them tucked in
            tick.set_y(-0.02)   # negative moves labels closer to axis; try -0.01 to -0.05

    for j in range(n, len(axes)):
        axes[j].axis("off")
    
    # reduce font size
    plt.rcParams.update({
        "font.size": 7,
        "axes.titlesize": 8,
        "axes.labelsize": 7,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "axes.linewidth": 0.7,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "legend.frameon": False,
        "svg.fonttype": "none",
    })

    fig.suptitle(
        "Participation-profile similarity (per-structure mean): within vs each other peptide (1 vs 5)",
        y=1.02,
        fontsize=10,
    )
    fig.savefig(out_p, dpi=600, bbox_inches="tight")
    fig.savefig(out_p.with_suffix('.svg'), dpi=600, bbox_inches="tight")
    plt.close(fig)


def _plot_participation_ecf(
    df_embedding: pd.DataFrame,
    ref_participation_array: list,  # kept for signature compatibility; ignored
    out_p: pathlib.Path,
    palette: dict,
    ncols: int = 3,
    show_points: bool = False,
    out_dir: pathlib.Path = None,
):
    """
    For each peptide (subplot), plot empirical CDF curves of Euclidean distances
    between that peptide's participation vectors and the participation vectors of
    OTHER peptides (grouped by other_peptide).

    Expected df_embedding columns:
      - peptide
      - participation_array (list/array-like)
      - other_peptide (optional; will be constructed if missing)

    Notes:
      - Reference is ignored entirely (ref_participation_array unused).
      - Distances are computed pairwise: (peptide == pep) vs (peptide == other).
      - If your df already has an 'other_peptide' column, it will be overwritten
        to enforce "other peptides" semantics.
    """
    if out_dir is None:
        out_dir = pathlib.Path(out_p).parent
    out_dir = pathlib.Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df_embedding = df_embedding.copy()

    required = {"peptide", "participation_array"}
    missing = required - set(df_embedding.columns)
    if missing:
        raise ValueError(f"df_embedding missing columns: {missing}")

    # --- Build a long-form table of distances: peptide vs other_peptide ---
    peptides = list(pd.unique(df_embedding["peptide"]))
    part_mat = np.vstack(df_embedding["participation_array"].to_numpy()).astype(float)

    # map peptide -> row indices in df_embedding
    idx_by_pep = {
        pep: np.flatnonzero(df_embedding["peptide"].to_numpy() == pep)
        for pep in peptides
    }

    rows = []
    for pep in peptides:
        idx_i = idx_by_pep.get(pep, np.array([], dtype=int))
        if idx_i.size == 0:
            continue
        A = part_mat[idx_i, :]  # (n_i, d)

        for other in peptides:
            if other == pep:
                continue  # ignore self
            idx_j = idx_by_pep.get(other, np.array([], dtype=int))
            if idx_j.size == 0:
                continue
            B = part_mat[idx_j, :]  # (n_j, d)

            # pairwise Euclidean distances between all A and all B
            # result shape (n_i, n_j)
            dists = np.linalg.norm(A[:, None, :] - B[None, :, :], axis=2).ravel()

            for d in dists:
                rows.append((pep, other, float(d)))

    plot_df = pd.DataFrame(rows, columns=["peptide", "other_peptide", "euclidean_distance"])
    if plot_df.empty:
        raise ValueError("No pairwise distances computed (check peptide groups / participation_array).")

    # --- Supplementary-figure style (clean, print-friendly) ---
    sns.set_theme(style="white", context="paper")
    plt.rcParams.update({
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "legend.frameon": False,
        "svg.fonttype": "none",
        "lines.solid_joinstyle": "round",
        "lines.solid_capstyle": "round",
    })

    n = len(peptides)
    ncols_eff = max(1, min(ncols, n))
    nrows = int(np.ceil(n / ncols_eff))

    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols_eff,
        figsize=(3.2 * ncols_eff, 2.6 * nrows),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    axes = np.atleast_1d(axes).ravel()

    # Global x limits (consistent across panels)
    x = plot_df["euclidean_distance"].to_numpy()
    x = x[np.isfinite(x)]
    x_min, x_max = np.min(x), np.max(x)
    pad = 0.03 * (x_max - x_min if x_max > x_min else 1.0)
    xlims = (x_min - pad, x_max + pad)

    def _ecdf_step(arr: np.ndarray):
        arr = np.asarray(arr, dtype=float)
        arr = arr[np.isfinite(arr)]
        if arr.size == 0:
            return np.array([]), np.array([])
        xs = np.sort(arr)
        ys = np.arange(1, xs.size + 1, dtype=float) / xs.size
        return xs, ys

    for i, pep in enumerate(peptides):
        ax = axes[i]
        sub = plot_df.loc[plot_df["peptide"] == pep].copy()

        # Sort curves by median distance (low -> high)
        med = (
            sub.groupby("other_peptide", as_index=False)["euclidean_distance"]
            .median()
            .sort_values("euclidean_distance", ascending=True)
        )
        order = med["other_peptide"].tolist()

        for other in order:
            arr = sub.loc[sub["other_peptide"] == other, "euclidean_distance"].to_numpy()
            xs, ys = _ecdf_step(arr)
            if xs.size == 0:
                continue

            color = palette.get(other, "#666666")

            ax.step(
                xs, ys,
                where="post",
                linewidth=1.25,
                color=color,
                alpha=0.95,
                label=other,
                zorder=2,
                solid_capstyle="round",
                solid_joinstyle="round",
            )

            if show_points:
                ax.plot(
                    xs, ys,
                    linestyle="None",
                    marker="o",
                    markersize=1.0,
                    alpha=0.18,
                    color=color,
                    zorder=1,
                )

        ax.set_title(
            pep,
            pad=6,
            bbox=dict(
                facecolor=palette.get(pep, "#ffffff"),
                edgecolor="black",
                linewidth=0.6,
                boxstyle="round,pad=0.22",
            ),
        )
        ax.set_xlim(*xlims)
        ax.set_ylim(0, 1.0)

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.tick_params(axis="both", which="major", direction="out")

        if i % ncols_eff == 0:
            ax.set_ylabel("Fraction ≤ distance")
        else:
            ax.set_ylabel("")

        if i >= (nrows - 1) * ncols_eff:
            ax.set_xlabel("Euclidean distance (participation profile)")
        else:
            ax.set_xlabel("")

        ax.legend(
            title=None,
            fontsize=6.5,
            loc="lower right",
            frameon=False,
            handlelength=2.4,
            handletextpad=0.6,
            borderaxespad=0.25,
        )

        ax.grid(
            True,
            axis="both",
            linewidth=1.0,
            alpha=0.25,
            color="gray",
            zorder=0,
            linestyle="--",
        )

    for j in range(n, len(axes)):
        axes[j].axis("off")

    fig.suptitle(
        "Empirical CDF of residue participation profile distances (vs other peptides)",
        y=1.02,
        fontsize=10,
    )

    # Save: respect out_p if provided; otherwise mirror your previous behavior
    save_path = out_p if out_p is not None else (out_dir / "residue_participation_profile_distances_ecdf.png")
    fig.savefig(save_path, dpi=600, bbox_inches="tight")
    plt.close(fig)


def main():
    ""
    # define paths
    script_dir = pathlib.Path(__file__).parent.resolve()
    df = pd.read_csv(script_dir / "run.csv")

    # GIP            250
    # Exendin-9      250
    # Semaglutide    250
    # GLP-1          250
    # Exendin-4      250
    # Tirzepatide    250
    peptides_to_keep = [
        "GIP",
        # "Exendin-9",
        "Semaglutide",
        "GLP-1",
        "Exendin-4",
        "Tirzepatide",
    ]
    df = df[df['peptide'].isin(peptides_to_keep)].copy()
    out_dir = script_dir / "peptide_dimerization_profile_similarities"
    out_dir.mkdir(exist_ok=True)

    # calc all vs all if not exists
    all_vs_all_distances_p = out_dir / "all_vs_all_peptide_dimerization_profile_euclidean_distances.csv"
    if not all_vs_all_distances_p.exists():
        calculate_peptide_dimerization_profile_distances(df, all_vs_all_distances_p)
    # load and set self-distances to nan
    euc_df = pd.read_csv(all_vs_all_distances_p, index_col=0)
    np.fill_diagonal(euc_df.values, np.nan)    

    # take the best (iptm) per peptide, seed, model (i.e. reduce on sample)
    df = df.loc[df.groupby(['peptide', 'seed', 'model'])['iptm'].idxmax()]
    # value counts
    print("Peptide counts after selecting best iptm per peptide, seed, model:")
    print(df['peptide'].value_counts())

    # TODO: consider running rmsd as well, takes 7h
    # # # calc all vs all rmsd if not exists
    # # all_vs_all_rmsd_distances_p = out_dir / "all_vs_all_peptide_dimerization_profile_rmsd_distances.csv"
    # # if not all_vs_all_rmsd_distances_p.exists():
    # #     calculate_peptide_dimerization_profile_rmsd_distances(df, all_vs_all_rmsd_distances_p)


    # run residue participation analysis
    plot_residue_participation_analysis_per_peptide(
        df=df,
        out_dir=out_dir,
        palette={**COLORS, "Exendin-9": "#888888"},
        ncols=1,
    )

    # for each peptide, get the distances to other peptides
    peptides = df['peptide'].unique()
    plot_df = []
    for peptide in peptides:
        for other_peptide in peptides:
            if peptide == other_peptide:
                continue
            plot_df.extend(_fill_in_plot_df(euc_df, df, peptide, other_peptide))
    plot_df = pd.DataFrame(plot_df)
    plot_df_p = out_dir / "peptide_dimerization_profile_euclidean_distances_for_plotting.csv"
    plot_df.to_csv(plot_df_p, index=False)

    # TODO: lastly, plot boxplot comparison, where we see dimerization euclidean-distance of each peptide to each other peptide
    print(f"Plotting dataframe saved to {plot_df_p}")
    plot_p = out_dir / "peptide_dimerization_profile_euclidean_distances_boxplot.png"
    plot_boxplot_dimerization_profile_distances(plot_df, plot_p, palette=COLORS)
    plot_sorted_median_boxplots_per_peptide(plot_df=plot_df,
                                            out_p=out_dir / "peptide_dimerization_profile_euclidean_distances_boxplots_per_peptide.png",
                                            palette=COLORS,
                                            ncols=1,
                                            )
    # set exendin-9 to gray instead of white for visibility
    plot_dimerization_cumalitive_euclidean_distance(
        plot_df=plot_df,
        out_p=out_dir / "peptide_dimerization_profile_euclidean_distances_cdf_per_peptide.png",
        palette={**COLORS, "Exendin-9": "#888888"},
        ncols=1,
        show_points=False,
    )





if __name__ == "__main__":
    """
    stem,seed,peptide,model,sample,iptm,id,pdb_p,dimer_distances,GIPR_pep_distances,GLP1R_pep_distances,receptor_bound
    unrelaxed_model_4_multimer_v3_pred_0,3,GIP,4,0,0.34351741635863414,seed3|GIP|unrelaxed_model_4_multimer_v3_pred_0,/Users/kcd635/HauserGroup Dropbox/luuk hd/GLP1R_GIPR_dimers/AF2-dimer-deepsampling/GLP1R_GIPR_dimers_originalseq/seed3/glp1r_human___gipr_human___GIP/unrelaxed_model_4_multimer_v3_pred_0.pdb,/Users/kcd635/HauserGroup Dropbox/luuk hd/GLP1R_GIPR_dimers/AF2-dimer-deepsampling/residue_distances/seed3|GIP|unrelaxed_model_4_multimer_v3_pred_0_B|C.csv,/Users/kcd635/HauserGroup Dropbox/luuk hd/GLP1R_GIPR_dimers/AF2-dimer-deepsampling/residue_distances/seed3|GIP|unrelaxed_model_4_multimer_v3_pred_0_C|D.csv,/Users/kcd635/HauserGroup Dropbox/luuk hd/GLP1R_GIPR_dimers/AF2-dimer-deepsampling/residue_distances/seed3|GIP|unrelaxed_model_4_multimer_v3_pred_0_B|D.csv,GLP1R
    unrelaxed_model_3_multimer_v3_pred_3,3,GIP,3,3,0.35198026346496414,seed3|GIP|unrelaxed_model_3_multimer_v3_pred_3,/Users/kcd635/HauserGroup Dropbox/luuk hd/GLP1R_GIPR_dimers/AF2-dimer-deepsampling/GLP1R_GIPR_dimers_originalseq/seed3/glp1r_human___gipr_human___GIP/unrelaxed_model_3_multimer_v3_pred_3.pdb,/Users/kcd635/HauserGroup Dropbox/luuk hd/GLP1R_GIPR_dimers/AF2-dimer-deepsampling/residue_distances/seed3|GIP|unrelaxed_model_3_multimer_v3_pred_3_B|C.csv,/Users/kcd635/HauserGroup Dropbox/luuk hd/GLP1R_GIPR_dimers/AF2-dimer-deepsampling/residue_distances/seed3|GIP|unrelaxed_model_3_multimer_v3_pred_3_C|D.csv,/Users/kcd635/HauserGroup Dropbox/luuk hd/GLP1R_GIPR_dimers/AF2-dimer-deepsampling/residue_distances/seed3|GIP|unrelaxed_model_3_multimer_v3_pred_3_B|D.csv,GIPR
    unrelaxed_model_3_multimer_v3_pred_2,3,GIP,3,2,0.38188280597861324,seed3|GIP|unrelaxed_model_3_multimer_v3_pred_2,/Users/kcd635/HauserGroup Dropbox/luuk hd/GLP1R_GIPR_dimers/AF2-dimer-deepsampling/GLP1R_GIPR_dimers_originalseq/seed3/glp1r_human___gipr_human___GIP/unrelaxed_model_3_multimer_v3_pred_2.pdb,/Users/kcd635/HauserGroup Dropbox/luuk hd/GLP1R_GIPR_dimers/AF2-dimer-deepsampling/residue_distances/seed3|GIP|unrelaxed_model_3_multimer_v3_pred_2_B|C.csv,/Users/kcd635/HauserGroup Dropbox/luuk hd/GLP1R_GIPR_dimers/AF2-dimer-deepsampling/residue_distances/seed3|GIP|unrelaxed_model_3_multimer_v3_pred_2_C|D.csv,/Users/kcd635/HauserGroup Dropbox/luuk hd/GLP1R_GIPR_dimers/AF2-dimer-deepsampling/residue_distances/seed3|GIP|unrelaxed_model_3_multimer_v3_pred_2_B|D.csv,GLP1R
    """
    main()
