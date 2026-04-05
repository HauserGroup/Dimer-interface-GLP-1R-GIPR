"""
"""
import pandas as pd 
import pathlib
import matplotlib.pyplot as plt
import seaborn as sns

script_dir = pathlib.Path(__file__).parent
df = pd.read_csv('code/orthologs/ncbi/MISTIC1_CMI/CMI_annotated.csv')
# only keep if 'TM' in 'segment'
df = df[df['segment'].str.contains('TM')]

# plot 'CMI' distribution per 'GPCR'
plot_p = script_dir / 'test.png'
plt.figure(figsize=(4, 6))
# plot distribution
sns.histplot(data=df, x='CMI', hue='GPCR', stat='density', common_norm=False, bins=20)
plt.title('CMI distribution per GPCR')
plt.xlabel('CMI')
plt.ylabel('Density')
plt.legend(title='GPCR', loc='upper right')
plt.tight_layout()
plt.savefig(plot_p, dpi=600)
plt.close()