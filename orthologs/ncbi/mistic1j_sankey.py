"""
Sankey plot of GLP1R-GIPR TM-to-TM co-evolution
------------------------------------------------
• Reads a CSV with columns …, segment1, segment2, MIvalue, …
• Aggregates MI per TM pair (mean by default; change to .sum() if desired)
• Builds a 14-node Sankey: 7 nodes for GLP1R (left), 7 for GIPR (right)
• Edge width ∝ MIvalue
"""

import pandas as pd
import plotly.graph_objects as go

# ------------------------------------------------------------------
# 1) LOAD YOUR DATA  (replace the filename with your real path)
# ------------------------------------------------------------------
df = pd.read_csv("code/orthologs/ncbi/MISTIC1_MI_pertubation/only_inter_tm_df.csv")
# keep only MI >6.5
df = df[df["MIvalue"] > 6.5]

# ------------------------------------------------------------------
# 2) AGGREGATE MI PER TM×TM COMBINATION
#    – mean MI emphasises "typical" strength
#    – swap to .sum() for total signal, or .max() for strongest contact
# ------------------------------------------------------------------
agg = (
    df.groupby(["segment1", "segment2"])["MIvalue"]
      .sum()                              # ← change here if needed
      .reset_index()
)

# ------------------------------------------------------------------
# 3) DEFINE NODES: GLP1R TM1-TM7 (indices 0-6)  →  GIPR TM1-TM7 (7-13)
# ------------------------------------------------------------------
tm_order = [f"TM{i}" for i in range(1, 8)]

left_nodes  = [f"GLP1R {tm}" for tm in tm_order]
right_nodes = [f"GIPR  {tm}" for tm in tm_order]   # double-space keeps labels aligned
nodes       = left_nodes + right_nodes

label_to_idx = {label: i for i, label in enumerate(nodes)}

# ------------------------------------------------------------------
# 4) BUILD LINK LISTS
# ------------------------------------------------------------------
sources, targets, values = [], [], []

for row in agg.itertuples(index=False):
    src_label = f"GLP1R {row.segment1}"
    tgt_label = f"GIPR  {row.segment2}"
    # Guard against pairs outside TM1-TM7 (optional)
    if src_label in label_to_idx and tgt_label in label_to_idx:
        sources.append(label_to_idx[src_label])
        targets.append(label_to_idx[tgt_label])
        values.append(row.MIvalue)

# ------------------------------------------------------------------
# 5) DRAW SANKEY
# ------------------------------------------------------------------
fig = go.Figure(
    data=go.Sankey(
        arrangement="snap",
        node=dict(
            pad=16,                # gap between nodes
            thickness=18,          # node bar width
            line=dict(width=0.5),
            label=nodes,
            hovertemplate="%{label}<extra></extra>",
        ),
        link=dict(
            source=sources,
            target=targets,
            value=values,
            hovertemplate=(
                "%{value:.3f} MI<br>"
                "<b>%{source.label}</b> → <b>%{target.label}</b>"
                "<extra></extra>"
            ),
        ),
    )
)

fig.update_layout(
    title="GLP1R – GIPR TM-to-TM mutual-information flow",
    font_size=12,
    margin=dict(l=20, r=20, t=40, b=20),
)

fig.show()
