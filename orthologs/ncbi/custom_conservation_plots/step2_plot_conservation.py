"""
1. From the jalview export, only take the residues of GLP1R.
"""
import pandas as pd
import pathlib as pl
import seaborn as sns
import matplotlib.pyplot as plt
sns.set_style("whitegrid")
sns.set_context('paper')




if __name__ == "__main__":
    script_dir = pl.Path(__file__).parent
    plot_p = script_dir / 'conservation_plots'
    plot_p.mkdir(exist_ok=True)


    dfs = ['GLP1R_jalview_annotation_fixed.csv',
           'GIPR_jalview_annotation_fixed.csv']
    for df in dfs:
        gpcr = df.split('_')[0]
        df_p = script_dir / df
        annot = pd.read_csv(df_p, index_col=0)
        for var in ['Conservation', 'Quality', 'Consensus', 'Occupancy', 'SHENKIN']:
            # plot quality of the alignment
            fig, ax = plt.subplots(figsize=(50, 5))
            sns.scatterplot(data=annot, x=annot.index, y=var, ax=ax)
            ax.set_xticks(annot.index)
            ax.set_xticklabels(annot['Residue'], rotation=0, fontsize=5)

            if min(annot[var]) < 0:
                raise ValueError(f"min {var} < 0")
            else:
                plt.ylim(0, None)

            # for tm4, do axvline
            for i, is_tm in zip(annot.index, annot['in_TM']):
                if is_tm:
                    ax.axvline(i, color='red', linestyle='--')

            # show interface
            for i, in_interface in zip(annot.index, annot['in_interface']):
                if in_interface:
                    ax.axvline(i, color='black', linestyle='--')
            
            plt.xlim(0, len(annot))
            plt.grid(axis='x')
            plt.savefig(plot_p / f"{gpcr}_{var}.png")
            plt.close()
