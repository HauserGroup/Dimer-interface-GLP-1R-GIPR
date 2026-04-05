"""
"""
import pathlib
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style='whitegrid')
sns.set_context('paper')

# load
script_dir = pathlib.Path(__file__).parent.absolute()
data_dir = script_dir / 'MISTIC1'
out_dir = script_dir / 'MISTIC1_MI'
out_dir.mkdir(exist_ok=True)

# sequence_number,amino_acid,protein_segment,display_generic_number
# therefore, we expect 929 residues in total
glp1r_info = pd.read_csv('code/glp1r_human_chain_info.csv')
glp1r_info['gpcr'] = 'GLP1R'
glp1r_info['gpcr_sequence_number'] = glp1r_info['sequence_number']
print(glp1r_info.shape) # (463, 4)
# chains
gipr_info = pd.read_csv('code/gipr_human_chain_info.csv')
gipr_info['gpcr'] = 'GIPR'
gipr_info['gpcr_sequence_number'] = gipr_info['sequence_number']
gipr_info['sequence_number'] = gipr_info['sequence_number'] + glp1r_info.shape[0]
assert max(gipr_info['sequence_number']) == 929, 'GIPR last residue sequence number is not 929'
assert len(set(gipr_info['sequence_number']).intersection(set(glp1r_info['sequence_number']))) == 0, 'GIPR and GLP1R sequence numbers overlap'
chain_info = pd.concat([glp1r_info, gipr_info], ignore_index=True)
chain_info.to_csv(out_dir / 'combined_MSA_chain_info.csv', index=False)
chain_info.set_index('sequence_number', inplace=True)

# load mutual information
mi_data = data_dir / 'MI_DATA'
mi = pd.read_csv(mi_data, sep='\t', comment='#', names=['resnum1', 'resname1', 'resnum2', 'resname2', 'MIvalue'])
mi['GPCR1'] = chain_info.loc[mi['resnum1'], 'gpcr'].values
mi['segment1'] = chain_info.loc[mi['resnum1'], 'protein_segment'].values
mi['GPCR2'] = chain_info.loc[mi['resnum2'], 'gpcr'].values
mi['segment2'] = chain_info.loc[mi['resnum2'], 'protein_segment'].values
mi['gpcr_sequence_number1'] = chain_info.loc[mi['resnum1'], 'gpcr_sequence_number'].values
mi['gpcr_sequence_number2'] = chain_info.loc[mi['resnum2'], 'gpcr_sequence_number'].values
print(min(mi['resnum1']), max(mi['resnum1'])) # 1-928
print(min(mi['resnum2']), max(mi['resnum2'])) # 2-929

# print top 5 MI values
print(mi.sort_values('MIvalue', ascending=False).head(5))
print(mi.sort_values('MIvalue', ascending=True).head())

# fix NaN problem
print(mi[mi['MIvalue'] == -999.99].shape) # (928, 9)
nan_values= mi[mi['MIvalue'] == -999.99]
nan_values.to_csv(out_dir / 'MI_NA_values.csv', index=False)
print(nan_values['resnum1'].unique())
print(nan_values['resnum2'].unique())
nan_value_count = nan_values.shape[0]
mi.loc[mi['MIvalue'] < -999, 'MIvalue'] = pd.NA

mi.to_csv(out_dir / 'MI_annotated.csv', index=False)


# plot MI distribution of all residues
plt.figure(figsize=(10, 5))
sns.histplot(mi['MIvalue'], bins=100, color='blue')
plt.xlabel('Mutual Information')
plt.ylabel('Frequency')
plt.title(f'MISTIC1 MI distribution (n={mi.shape[0]})')
plt.tight_layout()
plt.text(x=max(mi['MIvalue']) - int(max(mi['MIvalue']) * 0.1), 
         y=0, verticalalignment='bottom', horizontalalignment='right',
         fontsize=12, color='red', s=f'NA values: {nan_value_count}')
plt.savefig(out_dir / 'MI_distribution.png')
plt.close()

# get an array with GLP1R x GIPR MI values
mi_matrix = mi.pivot(index='resnum1', columns='resnum2', values='MIvalue')
mi_matrix.to_csv(out_dir / 'MI_matrix.csv')
print(mi_matrix.shape) # (928, 928)
print(max(mi_matrix.values.flatten())) # 33.22
print(min(mi_matrix.values.flatten())) # -999.99

# plot as a heatmap
plt.figure(figsize=(10, 5))
sns.heatmap(mi_matrix, cmap='viridis', cbar_kws={'label': 'MI', 'orientation': 'horizontal', 'shrink': 0.5},
            square=True)
plt.xlabel('Residue number')
plt.ylabel('Residue number')
plt.title('MISTIC1 MI')
plt.tight_layout()
plt.savefig(out_dir / 'MI_heatmap.png')