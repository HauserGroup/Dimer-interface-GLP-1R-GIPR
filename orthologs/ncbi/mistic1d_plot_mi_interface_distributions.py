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
out_dir = script_dir / 'MISTIC1_MI_interface2'
out_dir.mkdir(exist_ok=True)

# read the MI_annotated .csv
# resnum1,resname1,resnum2,resname2,MIvalue,GPCR1,segment1,GPCR2,segment2,gpcr_sequence_number1,gpcr_sequence_number2
mi = pd.read_csv(script_dir / 'MISTIC1_MI_interface/MI_inter_annotated.csv')

glp1r_segments = mi[mi['GPCR1'] == 'GLP1R']['segment1'].unique()
gipr_segments = mi[mi['GPCR2'] == 'GIPR']['segment2'].unique()
assert set(glp1r_segments) == set(gipr_segments), 'GLP1R and GIPR segments are not the same'

# get median MI value for each segment pair
mi['segment_pair'] = mi['GPCR1']+mi['segment1'] + 'x' + mi['GPCR2']+mi['segment2']

# only keep TMs
mi = mi[mi['segment1'].str.contains('TM') & mi['segment2'].str.contains('TM')]

# for each pair, plot the mean and std of the MI values
fig, ax = plt.subplots(figsize=(14, 7))
sns.violinplot(x='segment_pair', y='MIvalue', data=mi, ax=ax, linewidth=0.5)
ax.set_ylabel('MI value')
ax.set_xlabel('Segment pair')
ax.set_title('MISTIC1 MI interface distributions')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(out_dir / 'MISTIC1_MI_interface_distributions.png')
plt.close()
