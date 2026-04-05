#!/usr/bin/env python
"""
Circular Chord (Circos-style) of GPCR–TM links
------------------------------------------------
• Reads only_inter_tm_df.csv
• (Optionally) keeps top-5 % MI rows
• Aggregates mean MI for every TM pair
• Builds a 14 × 14 adjacency matrix
• Renders an interactive chord diagram

Requires:
    pip install pandas plotly
"""

import numpy as np
import pandas as pd
import plotly.figure_factory as ff
import pathlib

# ─── CONFIG ──────────────────────────────────────────────────────────────
SCRIPT_DIR = pathlib.Path(__file__).parent
CSV_PATH      = SCRIPT_DIR / 'MISTIC1_MI_pertubation' / 'only_inter_tm_df.csv'
KEEP_TOP_PCT  = 50         # 0‒100 ; set 0 to skip thresholding
OUT_HTML      = SCRIPT_DIR / 'figures' / 'TM_pairs_chord_diagram.html'

# 14 circular sectors in this exact order
TM_LEFT  = [f"GLP1R_TM{i}" for i in range(1, 8)]
TM_RIGHT = [f"GIPR_TM{i}"  for i in range(1, 8)]
CATEGORIES = TM_LEFT + TM_RIGHT                       # len = 14
INDEX = {tm: i for i, tm in enumerate(CATEGORIES)}    # name → position
# ─────────────────────────────────────────────────────────────────────────

def load_filtered(csv_path: pathlib.Path, top_pct: int) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    if top_pct > 0:
        thresh = df['MIvalue'].quantile(1 - top_pct / 100)
        df = df[df['MIvalue'] >= thresh]
        print(f"Keeping top {top_pct}% (MI ≥ {thresh:.3f}) → {len(df)} rows")
    return df

def build_matrix(df: pd.DataFrame) -> np.ndarray:
    """14 × 14 symmetric matrix of mean MI per TM pair."""
    # Aggregate by TM pair (order doesn’t matter - make it symmetric)
    agg = (
        df.groupby(['segment1', 'segment2'])['MIvalue']
          .mean()
          .reset_index()
    )

    mat = np.zeros((14, 14))
    for _, row in agg.iterrows():
        left  = f"GLP1R_{row.segment1}"
        right = f"GIPR_{row.segment2}"
        i, j  = INDEX[left], INDEX[right]
        mat[i, j] = mat[j, i] = row.MIvalue            # symmetry
    return mat

import numpy as np
import plotly.graph_objects as go

def plot_chord(matrix: np.ndarray, labels: list[str], out_html: str) -> None:
    """
    Draw a circular chord diagram with Plotly ≥ 6 (no figure_factory).
    `matrix` must be symmetric (N×N), `labels` length N.
    """
    N = len(labels)
    totals = matrix.sum(axis=1)
    gap = 2 * np.pi * 0.005           # 0.5 % gap between sectors
    full_circle = 2 * np.pi - N * gap
    angles = full_circle * totals.cumsum() / totals.sum()
    starts = np.concatenate(([0.0], angles[:-1])) + np.arange(N) * gap
    ends   = angles + np.arange(N) * gap

    # ----- node arcs -----------------------------------------------------
    node_traces = []
    for i, (lab, a0, a1) in enumerate(zip(labels, starts, ends)):
        th = np.linspace(a0, a1, 30)
        r  = np.ones_like(th)
        x, y = r * np.cos(th), r * np.sin(th)
        node_traces.append(
            go.Scatter(
                x=x, y=y, mode="lines",
                line=dict(width=20, color="lightblue"),
                hoverinfo="text",
                text=f"{lab}<br>Total MI: {totals[i]:.2f}",
                showlegend=False,
            )
        )

    # ----- ribbons (edges) ----------------------------------------------
    edge_shapes = []
    max_w = matrix.max()
    for i in range(N):
        for j in range(i + 1, N):
            w = matrix[i, j]
            if w == 0:
                continue
            # use the mid-angle of each node for the ribbon ends
            a_i = (starts[i] + ends[i]) / 2
            a_j = (starts[j] + ends[j]) / 2
            xi, yi = np.cos(a_i), np.sin(a_i)
            xj, yj = np.cos(a_j), np.sin(a_j)
            edge_shapes.append(
                dict(
                    type="path",
                    path=f"M {xi},{yi} "
                         f"Q 0,0 {xj},{yj}",
                    line=dict(
                        width=0.5 + 6 * w / max_w,
                        color="rgba(100,100,100,0.6)"
                    )
                )
            )

    # ----- assemble figure ----------------------------------------------
    fig = go.Figure(node_traces)
    fig.update_layout(
        title="Inter-TM mutual information (chord diagram)",
        width=800, height=800,
        shapes=edge_shapes,
        xaxis=dict(showgrid=False, zeroline=False, visible=False),
        yaxis=dict(showgrid=False, zeroline=False, visible=False),
        margin=dict(l=40, r=40, t=60, b=40)
    )
    fig.write_html(out_html)
    print(f"Chord diagram saved → {out_html}")


def main() -> None:
    df = load_filtered(CSV_PATH, KEEP_TOP_PCT)
    M  = build_matrix(df)
    print(f"Adjacency matrix shape: {M.shape}")
    print(f"Adjacency matrix sum: {M.sum():.2f}")
    print(f"Adjacency matrix max: {M.max():.2f}")
    plot_chord(M, CATEGORIES, OUT_HTML)

if __name__ == '__main__':
    main()
