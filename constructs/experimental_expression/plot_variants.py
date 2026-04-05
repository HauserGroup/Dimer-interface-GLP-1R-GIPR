"""
1. convert excel to csv
2. filter for non loss of function variants
3. plot the change in surface expression for each variant
"""
import pathlib
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def replace_newline_with_space(df):
    """
    # fixed names:
    Index(['rsID', 'Chromosome', 'Position', 'Reference', 'Alternate', 'Source',
        'Filters - exomes', 'Filters - genomes', 'Transcript Consequence',
        'Annotation', 'Flags', 'Allele Count', 'Allele Number',
        'Allele Frequency', 'Log Δ Surface expression', 'SE', 'N',
        'Log Δ Mini-Gs', 'SE.1', 'N.1', 'LoF vs WT',
        'LoF, controlled for surface expression',
        'Log Δ GLP-1R endocytosis response', 'SE.2', 'N.2', 'LoF vs WT.1',
        'LoF, controlled for surface expression.1'],
        dtype='object')
    """
    # for each column, replace newline with space
    col_names = df.columns
    new_names = [n.replace('\n', ' ') for n in col_names]
    new_names = [n.strip() for n in new_names]
    new_names = [n.replace('  ', ' ') for n in new_names]

    df.columns = new_names
    return df


def get_table_df():
    """ Save as csv, and return dataframe
    """
    # convert excel to csv
    script_dir = pathlib.Path(__file__).parent.absolute()
    table_p = script_dir / 'table12.xlsx'
    table_s = script_dir / 'table12.csv'
    if not table_s.exists():
        df = pd.read_excel(table_p, index_col=0)
        df = replace_newline_with_space(df)
        df.to_csv(table_s)
    else:
        df = pd.read_csv(table_s, index_col='Protein consequence')
    return df 


def plot_change_in_surface_expression_vs_mini_gs_recruitment(df, save_p):
    """
    """
    # plot change in surface expression vs mini-Gs recruitment
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111)
    sns.scatterplot(x='Log Δ Surface expression', y='Log Δ Mini-Gs', data=df
                    , ax=ax)
    ax.set_xlabel('Log Δ Surface expression')
    ax.set_ylabel('Log Δ Mini-Gs')
    ax.set_title('Change in surface expression vs Mini-Gs recruitment')
    
    # xlim -2 to 1
    # ylim -2 to 0.5
    ax.set_xlim(-2, 1)
    ax.set_ylim(-2, 0.5)
    plt.tight_layout()
    plt.savefig(save_p)
    plt.close()


def get_aa_pos(s:str):
    """
    """
    n = s[1:-1]
    # check numeric
    assert n.isnumeric(), f'Error: {s} is not numeric'
    return int(n)


def get_tm4_range():
    """ tm4 range of GLP-1R 

    s261 to e292
    """
    return list(range(261, 292+1))


def check_sequence_matches(df):
    """
    """
    df['aa_pos'] = df.index.map(get_aa_pos)

    # https://www.uniprot.org/uniprotkb/P43220/entry#sequences
    glp1r_seq = "MAGAPGPLRLALLLLGMVGRAGPRPQGATVSLWETVQKWREYRRQCQRSLTEDPPPATDLFCNRTFDEYACWPDGEPGSFVNVSCPWYLPWASSVPQGHVYRFCTAEGLWLQKDNSSLPWRDLSECEESKRGERSSPEEQLLFLYIIYTVGYALSFSALVIASAILLGFRHLHCTRNYIHLNLFASFILRALSVFIKDAALKWMYSTAAQQHQWDGLLSYQDSLSCRLVFLLMQYCVAANYYWLLVEGVYLYTLLAFSVLSEQWIFRLYVSIGWGVPLLFVVPWGIVKYLYEDEGCWTRNSNMNYWLIIRLPILFAIGVNFLIFVRVICIVVSKLKANLMCKTDIKCRLAKSTLTLIPLLGTHEVIFAFVMDEHARGTLRFIKLFTELSFTSFQGLMVAILYCFVNNEVQLEFRKSWERWRLEHLHIQRDSSMKPLKCPTSSLSSGATAGSSMYTATCQASCS"

    for i, row in df.iterrows():
        aa_pos = row['aa_pos']
        wt_aa_according_to_df = i[0]
        wt_aa = glp1r_seq[aa_pos-1]
        assert wt_aa_according_to_df == wt_aa, f'Error: {wt_aa_according_to_df} != {wt_aa} at position {aa_pos}'
    print("All amino acids match")


def plot_GLP1R_variants(df, save_p, xlim=None, label_threshold=-0.50, figsize=(10, 6)):
    """
    # color on LOF vs WT
    # label each point with the variant name (protein consequence)

    # plot change in surface expression on y-axis
    # on x-axis, put the amino_acid_position

    """
    # stat dict for LOF
    print(df['LoF vs WT'].value_counts())

    # get amino acid position (apply to index)
    df['aa_pos'] = df.index.map(get_aa_pos)

    # plot change in surface expression vs amino acid position
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    sns.scatterplot(x='aa_pos', y='Log Δ Surface expression', data=df
                    , ax=ax, hue='LoF vs WT', palette='viridis_r')
    ax.set_xlabel('Amino acid position')
    ax.set_ylabel('Log Δ Surface expression')
    ax.set_title('Change in surface expression for GLP1R variants')

    # add TM4 range (grey opaque)
    tm4_range = get_tm4_range()
    ax.axvspan(tm4_range[0], tm4_range[-1], color='grey', alpha=0.2)
    # add tm4 label
    ax.text(np.mean(tm4_range), -1.5, 'TM4', fontsize=12, ha='center')
    
    # label if change in surface expression is less than -0.50
    for i, row in df.iterrows():
        if row['Log Δ Surface expression'] < label_threshold:
            # if not in xlim, skip
            if xlim is not None:
                if row['aa_pos'] < xlim[0] or row['aa_pos'] > xlim[1]:
                    continue
            ax.text(row['aa_pos'], row['Log Δ Surface expression'], i, fontsize=6, ha='center',
                    va='bottom', rotation=25)

    # report variant stats for just TM4
    tm4_df = df[df['aa_pos'].isin(tm4_range)]
    print_cols = ['aa_pos', 'Log Δ Surface expression', 'LoF vs WT']
    print_df = tm4_df[print_cols]
    print(print_df)
    
    print("TM4 variant stats GLP1R", "="*20)
    print(tm4_df['LoF vs WT'].value_counts())

    if xlim is not None:
        ax.set_xlim(xlim)

    plt.tight_layout()
    plt.savefig(save_p)
    plt.close()


def report_on_useful_tm4_variants(df, surface_expression_threshold=-0.7):
    # get only useful variants (non LOF and delta log surface expression of -0.25 or less)
    # also, only TM4
    tm4_range = get_tm4_range()
    df['aa_pos'] = df.index.map(get_aa_pos)
    # filter for TM4
    tm4_df = df[df['aa_pos'].isin(tm4_range)]

    # meets surface expression threshold
    tm4_df = tm4_df[tm4_df['Log Δ Surface expression'] >= surface_expression_threshold]
    
    # sort on surface expression
    tm4_df = tm4_df.sort_values('Log Δ Surface expression', ascending=False)
    for i, row in tm4_df.iterrows():
        report = f"{i}\t{round(row['Log Δ Surface expression'], 3)}\t{row['LoF vs WT']}\t{round(row['Log Δ Mini-Gs'], 3)}"
        print(report)

    # color tm4 white
    tm4_range = get_tm4_range()
    for i in tm4_range:
        print(f'color white, resi {i} and chain A')
    # report in a pymol friendly format
    get_label_com = lambda s: f'label (chain A and resi {s[1:-1]} and name CA), "{s}"'
    get_color_com = lambda s: f'color gray40, resi {s} and chain A'
    get_color_com_red = lambda s: f'color red, resi {s} and chain A'
    get_sticks_com = lambda s: f'show sticks, resi {s} and chain A'

    for i, row in tm4_df.iterrows():
        aa_pos = i[1:-1]
        print(get_label_com(i))
        # check if surface expression is less than -0.5
        if row['Log Δ Surface expression'] < -0.5:
            print(get_color_com_red(aa_pos))
        else:
            print(get_color_com(aa_pos))

        print(get_sticks_com(aa_pos))
        # set opacity for residue at 20%


def plot_tm4_surface_expression_vs_mini_gs_recruitment(df, save_p):
    """
    """
    # get only TM4
    tm4_range = get_tm4_range()
    df['aa_pos'] = df.index.map(get_aa_pos)
    tm4_df = df[df['aa_pos'].isin(tm4_range)]

    # plot change in surface expression vs mini-Gs recruitment
    fig = plt.figure(figsize=(5,5))
    ax = fig.add_subplot(111)
    sns.scatterplot(x='Log Δ Surface expression', y='Log Δ Mini-Gs', data=tm4_df,
    # color on LOF vs WT
                hue='LoF vs WT', palette='viridis_r'
    )
    ax.set_xlabel('Log Δ Surface expression')
    ax.set_ylabel('Log Δ Mini-Gs')
    ax.set_title('Change in surface expression vs Mini-Gs recruitment\nfor TM4 variants')

    # add labels
    for i, row in tm4_df.iterrows():
        ax.text(row['Log Δ Surface expression'], row['Log Δ Mini-Gs'], 
                i, fontsize=6, ha='center',
                va='bottom', rotation=25)
        

    # save
    plt.tight_layout()
    plt.savefig(save_p)
    plt.close()


def run_main():
    """
    """
    script_dir = pathlib.Path(__file__).parent.absolute()
    plot_dir = script_dir / 'plots'
    plot_dir.mkdir(exist_ok=True)

    # convert excel to csv
    df = get_table_df()
    # check_sequence_matches(df.copy())

    # # we are using the Mini-Gs data, so drop endocytosis data
    # cols_to_drop = ['Log Δ GLP-1R endocytosis response', 'SE.2', 'N.2', 'LoF vs WT.1', 'LoF, controlled for surface expression.1']
    # df = df.drop(cols_to_drop, axis=1)

    # # check if correct data by remaking figure
    # plot_change_in_surface_expression_vs_mini_gs_recruitment(df.copy(),  plot_dir / 'change_in_surface_expression_vs_mini_gs_recruitment.png')

    # # for each non LOF variant, plot the change in surface expression
    # print("Variant stats GLP1R", "="*20)
    # plot_GLP1R_variants(df.copy(), plot_dir/'GLP1R_variants_labelled.svg')
    # plot_GLP1R_variants(df.copy(), plot_dir/'GLP1R_variants_labelled_only_tm4.svg', xlim=(260, 293), label_threshold=float('inf'),
    #                     figsize=(5, 5))

    # get only useful variants (non LOF and delta log surface expression of -0.25 or less)
    report_on_useful_tm4_variants(df.copy(), 
                                  surface_expression_threshold=-.5)

    # plot for tm4 the change in surface expression vs the change in mini-Gs recruitment
    plot_tm4_surface_expression_vs_mini_gs_recruitment(df.copy(), plot_dir / 'tm4_surface_expression_vs_mini_gs_recruitment.svg')






if __name__ == '__main__':
    run_main()

    # color white, resi 261-292 and chain A
    # set cartoon_transparency, 0.7, chain B
    # distance hydrophobic_interactions, (chain A and resi 261-292), (chain B ), 5
    # make labels bigger: 
    # set label_size, 30 
    # set cartoon_transparency, 0.0, chain B