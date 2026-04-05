"""
load only_inter_tm_df.csv

Take top 5% of MIvalue

resnum1,resname1,resnum2,resname2,MIvalue,GPCR1,segment1,GPCR2,segment2,gpcr_sequence_number1,gpcr_sequence_number2,Construct,ConstructImpact,IntraPair
164,A,620,L,22.064512,GLP1R,TM1,GIPR,TM1,164,157,,,False
192,L,782,G,21.271887,GLP1R,TM2,GIPR,TM5,192,319,,,False

Take unique pairs of 'GPCR1,segment1,GPCR2,segment2'
"""
"""
Analyse inter-TM mutual information.

• Load only_inter_tm_df.csv
• Keep the top 5 % of MIvalue
• Plot a 7 × 7 heat-map and a Sankey flow diagram
"""
import pathlib
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import networkx as nx


SEGMENTS = ["TM1", "TM2", "TM3", "TM4", "TM5", "TM6", "TM7"]


def load_and_filter(csv_path: pathlib.Path, quantile) -> pd.DataFrame:
    """Read CSV and retain the top-5 % MI rows."""
    df = pd.read_csv(csv_path)
    # keep unique pairs
    df = df.drop_duplicates(subset=["GPCR1", "segment1", "resnum1", "GPCR2", "segment2", "resnum2"])
    top5 = df["MIvalue"].quantile(quantile)
    print(f"Top-{quantile}% MIvalue threshold: {top5:0.3f}")
    return df[df["MIvalue"] >= top5]


def make_heatmap(df: pd.DataFrame, quantile, out: pathlib.Path, func=pd.DataFrame.mean) -> None:
    """7 × 7 heat-map: mean MI per TM–TM pair."""
    if func == pd.DataFrame.sum:
        pivot = (
            df.groupby(["segment1", "segment2"])["MIvalue"]
            .sum()          # Series (index: MultiIndex segment1/segment2)
            .unstack()       # → DataFrame, rows = segment1, cols = segment2
            .reindex(index=SEGMENTS, columns=SEGMENTS)  # enforce 7×7 order
            .fillna(0)       # optional: make missing pairs 0 (or NaN if you prefer)
        )
    
    elif func == pd.DataFrame.mean:
        pivot = (
            df.groupby(["segment1", "segment2"])["MIvalue"]
            .mean()          # Series (index: MultiIndex segment1/segment2)
            .unstack()       # → DataFrame, rows = segment1, cols = segment2
            .reindex(index=SEGMENTS, columns=SEGMENTS)  # enforce 7×7 order
            .fillna(0)       # optional: make missing pairs 0 (or NaN if you prefer)
        )
    else:
        raise ValueError(f"Unknown aggregation function: {func}")

    plt.figure(figsize=(7, 6))
    sns.heatmap(
        pivot,
        cmap="viridis",
        linewidths=0.5,
        square=True,
        cbar_kws={"label": "Mean MI"},
    )
    # add counts to heatmap
    for i, row in enumerate(pivot.iterrows()):
        for j, val in enumerate(row[1]):
            plt.text(
                j + 0.5,
                i + 0.5,
                f"{val:.2f}",
                ha="center",
                va="center",
                color="white" if val > 0 else "black",
            )
    plt.title(f"Inter-TM mutual information (top {quantile:.0%})\nnumber of pairs: {len(df)}")
    plt.ylabel("GPCR 1 segment")
    plt.xlabel("GPCR 2 segment")
    plt.tight_layout()
    plt.savefig(out, dpi=300)
    plt.savefig(out.with_suffix(".svg"))
    plt.close()
    print(f"Heat-map saved → {out}")


def make_sankey(df: pd.DataFrame, out_html: pathlib.Path) -> None:
    """Interactive Sankey diagram of TM–TM links."""
    # Aggregate to mean MI per unique TM pair
    agg = (
        df.groupby(["segment1", "segment2"])["MIvalue"]
        .mean()
        .reset_index()
    )

    # Sankey nodes
    nodes = SEGMENTS + [f"{s}_2" for s in SEGMENTS]  # suffix _2 for GPCR2 side
    node_indices = {name: i for i, name in enumerate(nodes)}

    # Links
    sources = [node_indices[row.segment1] for _, row in agg.iterrows()]
    targets = [node_indices[f"{row.segment2}_2"] for _, row in agg.iterrows()]
    values = agg["MIvalue"].tolist()
    labels = nodes

    sankey = go.Sankey(
        node=dict(label=labels, pad=15, thickness=15),
        link=dict(source=sources, target=targets, value=values),
    )

    fig = go.Figure(sankey)
    fig.update_layout(
        title_text="Inter-TM connectivity (top 5 %)",
        font_size=12,
        width=900,
        height=600,
    )
    fig.write_html(out_html)
    print(f"Sankey diagram saved → {out_html} (open in a browser)")



# ──────────────────────────────────────────────────────────────────────────
# 3) Bipartite force-directed network
# ──────────────────────────────────────────────────────────────────────────
def make_bipartite_graph(df: pd.DataFrame, out_html: pathlib.Path) -> None:
    """Interactive bipartite force-directed graph of TM-TM links."""
    agg = (
        df.groupby(["segment1", "segment2"])["MIvalue"]
          .mean()
          .reset_index()
    )

    # Build a bipartite graph
    G = nx.Graph()
    left_nodes  = [f"1_{s}" for s in SEGMENTS]   # GPCR1 side
    right_nodes = [f"2_{s}" for s in SEGMENTS]   # GPCR2 side
    G.add_nodes_from(left_nodes,  bipartite=0, label=[s[2:] for s in left_nodes])
    G.add_nodes_from(right_nodes, bipartite=1, label=[s[2:] for s in right_nodes])

    for _, row in agg.iterrows():
        G.add_edge(
            f"1_{row.segment1}",
            f"2_{row.segment2}",
            weight=row.MIvalue
        )

    # Force-directed layout (weighted spring)
    pos = nx.spring_layout(G, k=0.8, iterations=200, weight="weight", seed=42)

    # Edges
    edge_x, edge_y = [], []
    for u, v, attrs in G.edges(data=True):
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        edge_x += [x0, x1, None]
        edge_y += [y0, y1, None]

    edge_trace = go.Scatter(
        x=edge_x, y=edge_y, mode="lines",
        line=dict(width=2, color="gray"),
        hoverinfo="none"
    )

    # Nodes
    node_x = [pos[n][0] for n in G.nodes()]
    node_y = [pos[n][1] for n in G.nodes()]
    labels = [G.nodes[n]["label"] for n in G.nodes()]

    node_trace = go.Scatter(
        x=node_x, y=node_y, mode="markers+text",
        text=labels, textposition="bottom center",
        marker=dict(size=20, line=dict(width=2), color="lightblue"),
        hoverinfo="text"
    )

    fig = go.Figure([edge_trace, node_trace])
    fig.update_layout(
        title="Bipartite TM network (top 5 %)",
        showlegend=False,
        width=900, height=600,
        margin=dict(l=20, r=20, t=50, b=20)
    )
    fig.write_html(out_html)
    print(f"Bipartite network saved → {out_html}")


def create_cytoscape_csv(df, quantile, save_p):
    df['Source'] = df['GPCR1'] + '_' + df['segment1'] + '_' + df['resnum1'].astype(str) + df['resname1']
    df['Target'] = df['GPCR2'] + '_' + df['segment2'] + '_' + df['resnum2'].astype(str) + df['resname2']
    df[['Source', 'Target', 'MIvalue']].to_csv(save_p, index=False)
    print(f"Cytoscape CSV saved → {save_p}")


def main() -> None:
    script_dir = pathlib.Path(__file__).parent
    csv_path = script_dir / "MISTIC1_MI_pertubation" / "only_inter_tm_df.csv"

    figs_dir = script_dir / "figures"
    figs_dir.mkdir(exist_ok=True)

    for quantile in (0.50, 0.80, 0.90, 0.95, 0.99):
        print(f"Quantile: {quantile}")
        df_top = load_and_filter(csv_path, quantile)
        # Save cytoscape CSV
        create_cytoscape_csv(df_top, quantile, figs_dir / f"tm_cytoscape_{quantile}.csv")

        # plot
        make_heatmap(df_top, quantile, figs_dir / f"tm_heatmap_{quantile}.png",
                     func=pd.DataFrame.sum)
        make_sankey(df_top, figs_dir / f"tm_sankey_{quantile}.html")
        make_bipartite_graph(df_top, figs_dir / f"tm_bipartite_{quantile}.html")

    # plot heatmap for MI value > 6.5
    df = load_and_filter(csv_path, 0.0)
    # calculate which quantile 6.5 belongs in
    gt_6_5 = df[df["MIvalue"] > 6.5].shape[0] / df.shape[0]
    print(f"6.5 is in the {gt_6_5:.2%} quantile")
    df_top = df[df["MIvalue"] > 6.5]

    make_heatmap(df_top, gt_6_5, figs_dir / "tm_heatmap_6.5.png", func=pd.DataFrame.sum)
    make_sankey(df_top, figs_dir / "tm_sankey_6.5.html")
    


if __name__ == "__main__":
    main()
