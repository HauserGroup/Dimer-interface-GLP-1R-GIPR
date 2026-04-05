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
out_dir = script_dir / 'MISTIC1_CMI'
out_dir.mkdir(exist_ok=True)

# load combined MSA chain info
chain_info = pd.read_csv(script_dir / 'MISTIC1_MI/combined_MSA_chain_info.csv')
chain_info.set_index('sequence_number', inplace=True)

# load cumulative mutual information
cmi_data = data_dir / 'CMI_DATA'
cmi = pd.read_csv(cmi_data, sep='\t', comment='#', names=['resnum', 'CMI'])
# based on resnum, add which GPCR
cmi['GPCR'] = chain_info.loc[cmi['resnum'], 'gpcr'].values
cmi['segment'] = chain_info.loc[cmi['resnum'], 'protein_segment'].values
cmi.to_csv(out_dir / 'CMI_annotated.csv', index=False)

# plot CMI
plt.figure(figsize=(10, 5))
sns.scatterplot(x='resnum', y='CMI', data=cmi, label='CMI')
plt.xlabel('Residue number')
plt.ylabel('CMI')
plt.title('MISTIC1 CMI')
plt.legend()
plt.tight_layout()
plt.savefig(out_dir / 'CMI.png')

