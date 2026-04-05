"""
"""
import pathlib
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as patches
sns.set(style='whitegrid')
sns.set_context('paper')

# load
script_dir = pathlib.Path(__file__).parent.absolute()
data_dir = script_dir / 'MISTIC1'
out_dir = script_dir / 'MISTIC1_MI_interface'
out_dir.mkdir(exist_ok=True)

# read the MI_annotated .csv
# resnum1,resname1,resnum2,resname2,MIvalue,GPCR1,segment1,GPCR2,segment2,gpcr_sequence_number1,gpcr_sequence_number2
mi = pd.read_csv(script_dir / 'MISTIC1_MI/MI_annotated.csv')
# keep only pairs where GPCR1 != GPCR2
mi = mi[mi['GPCR1'] != mi['GPCR2']]
mi.to_csv(out_dir / 'MI_inter_annotated.csv', index=False)


# get segment coloring
# assert gipr segments are the same as glp1r segments
glp1r_segments = mi[mi['GPCR1'] == 'GLP1R']['segment1'].unique()
gipr_segments = mi[mi['GPCR2'] == 'GIPR']['segment2'].unique()
assert set(glp1r_segments) == set(gipr_segments), 'GLP1R and GIPR segments are not the same'
segments = glp1r_segments
segment_color_dict = dict()
for i, segment in enumerate(segments):
    if 'TM' in segment:
        segment_color_dict[segment] = 'darkblue'
    elif 'H8' in segment:
        segment_color_dict[segment] = 'lightblue'
    else:
        segment_color_dict[segment] = 'darkgreen'

# get indexes -> segments so we can color them later
mi_pivot = mi.pivot('gpcr_sequence_number1', 'gpcr_sequence_number2', 'MIvalue')
glp1r_indexes = mi_pivot.index
glp1r_index_to_segment = dict()
for glp1r_index in glp1r_indexes:
    segment = mi[mi['gpcr_sequence_number1'] == glp1r_index]['segment1'].values[0]
    glp1r_index_to_segment[glp1r_index] = segment
gipr_indexes = mi_pivot.columns
gipr_index_to_segment = dict()
for gipr_index in gipr_indexes:
    segment = mi[mi['gpcr_sequence_number2'] == gipr_index]['segment2'].values[0]
    gipr_index_to_segment[gipr_index] = segment
glp1r_colors = [segment_color_dict[glp1r_index_to_segment[glp1r_index]] for glp1r_index in glp1r_indexes]
gipr_colors = [segment_color_dict[gipr_index_to_segment[gipr_index]] for gipr_index in gipr_indexes]

# plot the heatmap (use clustermap to add the colors, but do not cluster)
fig, ax = plt.subplots(figsize=(7, 7))
cmap = sns.clustermap(mi_pivot,
               row_colors=gipr_colors, col_colors=gipr_colors, 
               row_cluster=False, col_cluster=False,
               cmap='coolwarm',
               # set label
               cbar_kws={'label': 'MI', 'orientation': 'horizontal'},
               xticklabels=False, yticklabels=False,
               colors_ratio=0.02,
)
# update clustermap to label x and y axis
cmap.ax_heatmap.set_xlabel('GLP1R', fontsize=12)
cmap.ax_heatmap.set_ylabel('GIPR', fontsize=12)
cmap.ax_heatmap.set_title('MI interface between GLP1R and GIPR', fontsize=14)
# update cbar
cbar = cmap.ax_heatmap.collections[0].colorbar
cbar.set_label('MI', fontsize=12)
cbar.ax.set_position([0.1, 0.9, 0.2, 0.025]) #x, y, width, height
plt.savefig(out_dir / 'MI_interface_heatmap.svg')
# disable title for png
cmap.ax_heatmap.set_title('')
cmap.ax_heatmap.set_xlabel('')
cmap.ax_heatmap.set_ylabel('')
plt.savefig(out_dir / 'MI_interface_heatmap.png', dpi=900)
plt.close()
