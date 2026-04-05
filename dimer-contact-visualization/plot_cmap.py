"""from a .pdb file, plot the cmap of the GLP1R x GIPR.

Plot the C-A distances per residue. Plot the segments
"""

import pathlib
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import PDB
import torch
import pandas as pd
import numpy as np

from settings import get_chains, get_protein_lens


def get_chain_info(chain_info_p):
    """reads the chain info from a csv file
    returns: pd.DataFrame
    with headers: 'sequence_number', 'amino_acid', 'protein_segment', 'display_generic_number'
    """
    df = pd.read_csv(chain_info_p)
    return df


def load_pdb_struct(pdb_p):
    """load a pdb file and return the structure"""
    parser = PDB.PDBParser()
    # struct name = filename without extension
    struct_name = pdb_p.stem
    structure = parser.get_structure(struct_name, pdb_p)
    return structure


def get_CA_atoms_per_chain(struct):
    """ """
    chains = struct.get_chains()
    chain_CA_atoms = {}
    for chain in chains:
        chain_CA_atoms[chain.id] = []
        for residue in chain:
            if residue.has_id("CA"):
                chain_CA_atoms[chain.id].append(residue["CA"])
    return chain_CA_atoms


def get_CA_coordinates_per_chain(struct):
    """ """
    chains = struct.get_chains()
    chain_CA_coord = {}
    for chain in chains:
        chain_CA_coord[chain.id] = torch.tensor(
            np.array(
                [
                    atom.get_coord()
                    for atom in chain.get_atoms()
                    if atom.get_id() == "CA"
                ]
            )
        )
    return chain_CA_coord


def plot_2D_CA_cmap(struct, chain_info, save_path, glp1r_chain="A", gipr_chain="B"):
    """Plot the CA distances between GLP1R and GIPR"""
    cmap = sns.color_palette("coolwarm", as_cmap=True)
    fig, ax = plt.subplots(figsize=(10, 10))

    chain_CA_atoms = get_CA_atoms_per_chain(struct)
    chain_CA_coord = {}
    for chain in chain_CA_atoms:
        chain_CA_coord[chain] = torch.tensor(
            [atom.get_coord() for atom in chain_CA_atoms[chain]]
        )
    chain_CA_cdist = torch.cdist(
        chain_CA_coord[glp1r_chain], chain_CA_coord[gipr_chain]
    )
    # plot the cmap
    sns.heatmap(
        chain_CA_cdist,
        cmap=cmap,
        ax=ax,
        cbar_kws={"label": "CA distance (Angstroms)"},
    )
    ax.set_xlabel(f"{gipr_chain} sequence number")
    ax.set_ylabel(f"{glp1r_chain} sequence number")
    ax.set_title("CA distances between GLP1R and GIPR")
    plt.savefig(save_path)
    plt.close()


def add_segments_to_x_axis(ax, chain_info):
    """Add the segments to the x-axis"""
    # keep track of segment
    old_segment_name = None
    segment_count = 0
    # for each AA
    for i, row in chain_info.iterrows():
        segment_name = row["protein_segment"]
        segment_index = row["sequence_number"]
        if old_segment_name is None:
            old_segment_name = segment_name
            old_segment_index = segment_index
        # if new segment
        if segment_name != old_segment_name:
            # plot the OLD segment
            color = ["grey" if segment_count % 2 == 0 else "white"][0]
            ax.axvspan(
                old_segment_index,
                segment_index,
                color=color,
                alpha=0.5,
            )
            ax.text(
                (segment_index + old_segment_index) / 2,
                0,
                old_segment_name,
                fontsize=8,
                rotation=90,
            )
            segment_count += 1
            old_segment_index = segment_index
            old_segment_name = segment_name

    # highlight the last segment
    ax.axvspan(
        old_segment_index,
        chain_info["sequence_number"].max(),
        color=["grey" if segment_count % 2 == 0 else "white"][0],
        alpha=0.5,
    )
    ax.text(
        (segment_index + old_segment_index) / 2,
        0,
        segment_name,
        fontsize=8,
        rotation=90,
    )


class HeavyAtomDistances:
    def __init__(self, structure, minimize_index, index_to_use) -> None:
        self.structure = structure
        self.minimize_index = minimize_index
        self.chain_to_minimize = list(structure.get_chains())[minimize_index]
        self.chain_to_use = list(structure.get_chains())[index_to_use]
        self.chains = structure.get_chains()

        self.atom_to_remove = "H"

    def group_atoms_per_residue(self, chain):
        """group atoms per residue"""
        chain_atoms = {}
        for residue in chain:
            chain_atoms[residue.id[1]] = []
            for atom in residue:
                if atom.get_id() != self.atom_to_remove:
                    chain_atoms[residue.id[1]].append(atom)
        return chain_atoms

    def get_distance_to_minimize_index(self):
        """get the distance to minimize index"""
        # get all the coordinates for the minimize index
        reference_coords = list()
        for residue in self.chain_to_minimize:
            for atom in residue:
                if atom.get_id() != self.atom_to_remove:
                    reference_coords.append(atom.get_coord())
        reference_coords = torch.tensor(np.array(reference_coords))

        # get the distances
        shortest_heavy_atom_distances = list()
        for residue in self.chain_to_use:
            residue_coords = list()
            for atom in residue:
                if atom.get_id() != self.atom_to_remove:
                    residue_coords.append(atom.get_coord())
            # calculate the distances to the reference
            residue_coords = torch.tensor(np.array(residue_coords))
            distances = torch.cdist(residue_coords, reference_coords)
            # get the smallest distance
            shortest_heavy_atom_distances.append(distances.min())

        return torch.tensor(shortest_heavy_atom_distances)


def plot_shortest_distances(
    struct_list,
    glp1r_chain_info,
    gipr_chain_info,
    save_path,
    glp1r_chain,
    gipr_chain,
    minimize_index,
    y_max=0,
    label_angst_thresh=0,
    xlim=None,
    distance_mode="CA",
):
    """
    GIPR is longer, so we will use its length
    x-axis = GLP1R residues
    y-axis = shortest distance between GLP1R and GIPR
    minimize_index: int,
                    0 to minimize GLP1R
                    1 to minimize GIPR
    distance_mode: str,
                    "CA" or "Heavy atoms"
    """
    # create empty fig
    fig, ax = plt.subplots(figsize=(10, 10))
    if minimize_index == 0:
        chosen_info = gipr_chain_info
        xlabel = f"GIPR sequence number"
    elif minimize_index == 1:
        chosen_info = glp1r_chain_info
        xlabel = f"GLP1R sequence number"

    # iterate over the structures and chain info
    for struct in struct_list:
        # get CA atoms
        if distance_mode == "CA":
            chain_CA_atoms = get_CA_atoms_per_chain(struct)
            # get CA coordinates per chain
            chain_CA_coord = get_CA_coordinates_per_chain(struct)
            # get the shortest distances
            chain_dist = torch.cdist(
                chain_CA_coord[glp1r_chain], chain_CA_coord[gipr_chain]
            )
            # minimize the distance the tensor on the gipr chain
            chain_dist = chain_dist.min(minimize_index).values
        elif distance_mode == "Heavy atoms":
            # get heavy atoms distances
            index_to_use = 1 if minimize_index == 0 else 0
            calculator = HeavyAtomDistances(struct, minimize_index, index_to_use)
            chain_dist = calculator.get_distance_to_minimize_index()
        # plot the shortest distances
        sns.scatterplot(
            x=chosen_info["sequence_number"],
            y=chain_dist,
            label=struct.get_id(),
            ax=ax,
        )
        if label_angst_thresh > 0:
            # check which points have yvals <= label_angst_thresh
            for i, y in enumerate(chain_dist):
                if y <= label_angst_thresh:
                    ax.text(
                        chosen_info["sequence_number"][i],
                        y,
                        f"{chosen_info['amino_acid'][i]}{chosen_info['sequence_number'][i]}  {chosen_info['display_generic_number'][i]}",
                        fontsize=6,
                        # move so that the text is not on the point
                        rotation=270,
                    )

    # set the labels and title
    ax.set_xlabel(xlabel)
    ax.set_ylabel("CA distance (Angstroms)")
    ax.set_title(f"Shortest {distance_mode} distances between GLP1R and GIPR")
    # put the legend on lower left, with a little bit of padding
    ax.legend(
        loc="lower left",
        bbox_to_anchor=(0.0, -0.1),
    )
    plt.tight_layout()
    # to the x-axis, add the GPCR segments
    add_segments_to_x_axis(ax, chosen_info)
    if y_max > 0:
        ax.set_ylim(0, y_max)
    # set xlim so that the x-axis starts at 0
    ax.set_xlim(0, gipr_chain_info["sequence_number"].max())
    if xlim is not None:
        ax.set_xlim(xlim)
    plt.savefig(save_path)
    plt.close()


if __name__ == "__main__":
    """
    # label the very frequent residues (opiod receptor residues, might be important but could be different)
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10497915/
    # very frequent (8–12x), moderate (3–7x), and rare (1–2x). Nineteen frequent residues were found: 110Met (ICL2), 160Ala (ECL2), Tyr3.48, Tyr3.51, Thr3.52, Ala3.55, Arg4.40, Arg4.41, Val4.44, Ile4.48, Leu4.52, Thr4.55, Pro4.59, Tyr5.41, Val5.45, Val5.49, Ile5.52, Val5.53, and Leu5.56. More detailed analyses of the various interfaces are 
    """
    GLP1_CHAIN_INF_P = "glp1r_human_chain_info.csv"
    GLP1_CHAIN_INF = get_chain_info(GLP1_CHAIN_INF_P)

    GIP_CHAIN_INF_P = "gipr_human_chain_info.csv"
    GIP_CHAIN_INF = get_chain_info(GIP_CHAIN_INF_P)

    GLP1R_GIPR_STRUCT = load_pdb_struct(
        pathlib.Path("data/AF2_PDBs/glp1r_human___gipr_human.pdb")
    )
    GLP1R_GIPR_GLP1_STRUCT = load_pdb_struct(
        pathlib.Path("data/AF2_PDBs/glp1r_human___gipr_human___glp1.pdb")
    )
    GLP1R_GIPR_GIP_STRUCT = load_pdb_struct(
        pathlib.Path("data/AF2_PDBs/glp1r_human___gipr_human___gip.pdb")
    )
    GLP1R_CHAIN = get_chains()["GLP1R"]
    GIPR_CHAIN = get_chains()["GIPR"]
    GLP1R_LEN = get_protein_lens()["GLP1R"]
    GIPR_LEN = get_protein_lens()["GIPR"]
    PLOT_DIR = pathlib.Path("plots/cmap")
    PLOT_DIR.mkdir(parents=True, exist_ok=True)

    # GIPR
    plot_shortest_distances(
        [GLP1R_GIPR_STRUCT],  # GLP1R_GIPR_GIP_STRUCT],
        GLP1_CHAIN_INF,
        GIP_CHAIN_INF,
        PLOT_DIR / "shortest_distances_minimize_GLP1R_HeavyAtoms.png",
        GLP1R_CHAIN,
        GIPR_CHAIN,
        minimize_index=0,
        distance_mode="Heavy atoms",
    )
    plot_shortest_distances(
        [GLP1R_GIPR_STRUCT],  # GLP1R_GIPR_GIP_STRUCT],
        GLP1_CHAIN_INF,
        GIP_CHAIN_INF,
        PLOT_DIR / "shortest_distances_minimize_GLP1R_HeavyAtoms_zoom.png",
        GLP1R_CHAIN,
        GIPR_CHAIN,
        y_max=15,
        minimize_index=0,
        label_angst_thresh=12,
        xlim=(120, 210),
        distance_mode="Heavy atoms",
    )
    plot_shortest_distances(
        [GLP1R_GIPR_STRUCT],  # GLP1R_GIPR_GIP_STRUCT],
        GLP1_CHAIN_INF,
        GIP_CHAIN_INF,
        PLOT_DIR / "shortest_distances_minimize_GLP1R_CA_zoom.png",
        GLP1R_CHAIN,
        GIPR_CHAIN,
        y_max=15,
        minimize_index=0,
        label_angst_thresh=12,
        xlim=(120, 210),
        distance_mode="CA",
    )

    # GLP1R
    plot_shortest_distances(
        [GLP1R_GIPR_STRUCT],  # GLP1R_GIPR_GIP_STRUCT],
        GLP1_CHAIN_INF,
        GIP_CHAIN_INF,
        PLOT_DIR / "shortest_distances_minimize_GIPR_HeavyAtoms.png",
        GLP1R_CHAIN,
        GIPR_CHAIN,
        minimize_index=1,
        distance_mode="Heavy atoms",
    )
    plot_shortest_distances(
        [GLP1R_GIPR_STRUCT],  # GLP1R_GIPR_GIP_STRUCT],
        GLP1_CHAIN_INF,
        GIP_CHAIN_INF,
        PLOT_DIR / "shortest_distances_minimize_GIPR_CA_zoom.png",
        GLP1R_CHAIN,
        GIPR_CHAIN,
        y_max=15,
        minimize_index=1,
        label_angst_thresh=12,
        xlim=(250, 320),
        distance_mode="CA",
    )
    plot_shortest_distances(
        [GLP1R_GIPR_STRUCT],  # GLP1R_GIPR_GIP_STRUCT],
        GLP1_CHAIN_INF,
        GIP_CHAIN_INF,
        PLOT_DIR / "shortest_distances_minimize_GIPR_HeavyAtoms_zoom.png",
        GLP1R_CHAIN,
        GIPR_CHAIN,
        y_max=15,
        minimize_index=1,
        label_angst_thresh=12,
        xlim=(250, 320),
        distance_mode="Heavy atoms",
    )


plot_shortest_distances(
    [
        GLP1R_GIPR_STRUCT,
        # GLP1R_GIPR_GIP_STRUCT,
        GLP1R_GIPR_GLP1_STRUCT,
    ],
    GLP1_CHAIN_INF,
    GIP_CHAIN_INF,
    PLOT_DIR / "test.png",
    GLP1R_CHAIN,
    GIPR_CHAIN,
    y_max=0,
    minimize_index=1,
    label_angst_thresh=0,
    xlim=None,
    distance_mode="CA",
)
