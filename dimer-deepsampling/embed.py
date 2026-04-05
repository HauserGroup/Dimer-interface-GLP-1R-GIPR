"""
Create embeddings (from the angstrom distances) to be used in a UMAP.

- use euclidean distance
- focus on short distance, clip to 0 after 8 Å
"""
import tqdm
import pandas as pd 
import pathlib
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from sklearn.linear_model import LinearRegression
import umap
import csv
sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.0)


def plot_embeddings(df, save_p, mode="tsne", dot_size=5, correct_for_len=True):
    """ df should contain feature_* columns and a 'label' column and a 'seq_id' column

    mode: "tsne" or "umap"
    """
    df = df.copy()

    # summarize features into umap
    features = df.filter(like='feature_').values

    if correct_for_len:
        length = df['sequence'].apply(len).values.reshape(-1, 1)
        # Regress out length from each feature
        corrected_features = []
        for i in range(features.shape[1]):
            model = LinearRegression().fit(length, features[:, i])
            residuals = features[:, i] - model.predict(length)
            corrected_features.append(residuals)
        corrected_features = np.array(corrected_features).T  # shape: (n_samples, n_features)
    else:
        corrected_features = features


    # plot the difference between the original features and the corrected features
    plt.figure(figsize=(3, 3))
    plt.scatter(features[:, 0], features[:, 1], s=1, alpha=0.5, label='Original Features')
    plt.scatter(corrected_features[:, 0], corrected_features[:, 1], s=1, alpha=0.5, label='Corrected Features')
    plt.xlabel('Feature 0')
    plt.ylabel('Feature 1')
    plt.title('Feature Correction by Length')
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_p.with_suffix('.feature_correction.png'), dpi=600)
    plt.close()

    if mode == "umap":
        reducer = umap.UMAP(random_state=42)
        embeddings = reducer.fit_transform(corrected_features)
    elif mode == "tsne":
        # to TSNE (dimensionality reduction)
        tsne = TSNE(n_components=2, random_state=42)
        embeddings = tsne.fit_transform(corrected_features)
    else:
        raise ValueError(f"Unknown mode: {mode}")

    # create a dataframe with the embeddings and labels
    embeddings_df = pd.DataFrame(embeddings, columns=['x', 'y'])
    embeddings_df['label'] = df['label'].values
    if 'seq_id' in df.columns:
        embeddings_df['id'] = df['seq_id'].values
    elif 'identifier' in df.columns:
        embeddings_df['id'] = df['identifier'].values
    else:
        raise ValueError("DataFrame must contain 'seq_id' or 'identifier' column for IDs.")

    # sort label, with nan last (if nan is present)
    label_order = embeddings_df['label'].dropna().unique().tolist()
    embeddings_df['label'] = embeddings_df['label'].fillna(pd.NA)
    if embeddings_df['label'].isna().any():
        label_order.append(pd.NA)


    # plot the embeddings
    fig, ax = plt.subplots(figsize=(3.0, 3.0))

    if len(label_order) > 8:
        print(f'Warning too many unique labels {len(label_order)}, colouring last gray')

    color_palette = sns.color_palette("colorblind", n_colors=8) # colorblind, deep, muted, cubehelix
    color_palette_dict = {}
    for i, label in enumerate(label_order):
        if i <= 7:
            color_palette_dict[label] = color_palette[i]
        else:
            color_palette_dict[label] = (0.8, 0.8, 0.8)

    embeddings_df['color'] = embeddings_df['label'].map(color_palette_dict)
    scatter = ax.scatter(
        embeddings_df['x'],
        embeddings_df['y'],
        c=embeddings_df['color'].tolist(),
        s=dot_size,
        alpha=0.7,
        edgecolors='none',
        rasterized=True, # rasterized means the points are rasterized but the axes are vectorized
    )

    # legend for each color
    label_value_counts = embeddings_df['label'].value_counts(dropna=False)
    handles = []
    for label, color in color_palette_dict.items():
        label_count = label_value_counts.get(label, 0)
        if pd.isna(label):
            display_label = 'NaN'
        else:
            display_label = str(label)
        display_label = f'{display_label} (n={label_count:,})'
        handles.append(plt.Line2D([0], [0], marker='o', color='w', label=display_label, markerfacecolor=color, markersize=5))
    plt.legend(handles=handles, loc='best', fontsize='small', markerscale=2, frameon=False)


    plt.title(f'{mode.upper()} of embeddings')
    plt.tight_layout()
    # set aspect ratio to 1:1
    plt.gca().set_aspect('equal', adjustable='box')
    # add xlabel, ylabel
    if mode == "umap":
        plt.xlabel('UMAP 1')
        plt.ylabel('UMAP 2')
    elif mode == "tsne":
        plt.xlabel('t-SNE 1')
        plt.ylabel('t-SNE 2')
    else:
        raise ValueError(f"Unknown mode: {mode}")
    
    # place the xtick and ytick close to the axes
    plt.tick_params(axis='both', which='both', length=0)  # remove ticks
    plt.subplots_adjust(left=0.05, right=0.95, top=0.98, bottom=0.05)
    
    plt.savefig(save_p, dpi=600)
    plt.savefig(save_p.with_suffix('.svg'))
    plt.close()

    # also save the embeddings with labels (tnse or umap coordinates)
    out_p = save_p.with_suffix(f'.{mode}_embeddings.csv')
    embeddings_df.to_csv(out_p, index=False)

# ---- plot embeddings end ----



def embed(df, out_p, angstrom_cutoff=12.0):
    """
    """
    expected_shape = (418, 410)
    with open(out_p, "w", newline="") as f:
        writer = None

        for row in tqdm.tqdm(df.itertuples(index=False), total=len(df), desc="Creating embeddings"):
            row_id = row.id
            row_label = row.peptide

            # load distance map
            dimer_p = pathlib.Path(row.dimer_distances)
            dimer_dist_maps = pd.read_csv(dimer_p, index_col=0)

            if dimer_dist_maps.shape != expected_shape:
                raise ValueError(
                    f"Unexpected shape for dimer distance map: {dimer_dist_maps.shape}, "
                    f"expected: {expected_shape}, file: {dimer_p}"
                )

            # flatten + clip
            dimer_dist_values = np.clip(
                dimer_dist_maps.values.ravel(),
                0,
                angstrom_cutoff,
            )

            # we wish to prioritize short distances, instead of raw distances, we use (angstrom_cutoff - distance)**2
            dimer_dist_values = (angstrom_cutoff - dimer_dist_values) ** 2

            # write header once
            if writer is None:
                header = (
                    ["seq_id", "label"]
                    + [f"feature_{i}" for i in range(len(dimer_dist_values))]
                )
                writer = csv.writer(f)
                writer.writerow(header)

            # write one row
            writer.writerow([row_id, row_label, *dimer_dist_values])
            

def main():
    script_dir = pathlib.Path(__file__).parent.resolve()
    df = pd.read_csv(script_dir / "run.scores.csv")
    out_p = script_dir / "run.embeddings.csv"
    
    embed(df, out_p, angstrom_cutoff=8.0)

    print('Reading embeddings for plotting...')
    embeddings = pd.read_csv(out_p)

    plot_dir = script_dir / "plots_embeddings"
    plot_dir.mkdir(exist_ok=True)

    print('Plotting umap...')
    plot_embeddings(
        embeddings,
        save_p = plot_dir / "embeddings_umap.png",
        mode="umap",
        dot_size=5,
        correct_for_len=False,
    )
    print("UMAP plot saved.")

    print('Plotting tsne...')
    plot_embeddings(
        embeddings,
        save_p = plot_dir / "embeddings_tsne.png",
        mode="tsne",
        dot_size=5,
        correct_for_len=False,
    )
    print("t-SNE plot saved.")
    

if __name__ == "__main__":
    main()

