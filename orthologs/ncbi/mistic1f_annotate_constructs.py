import pathlib
import pandas as pd


def gipr_constructs():
    """
    # constructs
    # 1: 136 143 147 151
    sele gipr_construct_1, resi 136+143+147+151 and gipr
    deselect
    # 2: 138 142 146 191
    sele gipr_construct_2, resi 138+142+146+191 and gipr
    deselect
    # 3: 136 138 142 143 146 147 151 191
    sele gipr_construct_3, resi 136+138+142+143+146+147+151+191 and gipr
    deselect
    """
    # construct number : [residues, ]
    return {
        5 : [136, 143, 147, 151],
        6 : [138, 142, 146, 191],
        7 : [136, 138, 142, 143, 146, 147, 151, 191],
    }


def glp1r_constructs():
    """
    # 1: 279 283 291
    sele glp1r_construct_1, resi 279+283+291 and glp1r
    deselect
    # 2: 283 286 290
    sele glp1r_construct_2, resi 283+286+290 and glp1r
    deselect
    # 3: 275 279 282
    sele glp1r_construct_3, resi 275+279+282 and glp1r
    deselect
    # 4: 275 279 282 283 286 290 291
    sele glp1r_construct_4, resi 275+279+282+283+286+290+291 and glp1r
    deselect
    """
    return {
        1 : [279, 283, 291],
        2 : [283, 286, 290],
        3 : [275, 279, 282],
        4 : [275, 279, 282, 283, 286, 290, 291],
    }


def get_construct_impact():
    return {
        'gipr' : {
            5:"None",
            6:"Strong",
            7:"Strong",
        },
        'glp1r': {
            1:"Medium",
            2:"Medium",
            3:"Strong",
            4:"Strong",
        }
    }


# paths
script_dir = pathlib.Path(__file__).parent
out_dir = script_dir / 'MISTIC1_MI_pertubation'
out_dir.mkdir(exist_ok=True)
df_p = 'code/orthologs/ncbi/MISTIC1_MI_interface/MI_inter_annotated.csv'
df = pd.read_csv(df_p)

# resnum1,resname1,resnum2,resname2,MIvalue,GPCR1,segment1,GPCR2,segment2,gpcr_sequence_number1,gpcr_sequence_number2,interface
# 1,M,2,A,12.502455,GLP1R,N-term,GLP1R,N-term,1,2,
# annotate construct numbers based on df.
df['Construct'] = pd.NA
df['ConstructImpact'] = pd.NA

# for each pair (intra or inter)
for i, row in df.iterrows():
    # for gpcr 1 or gpcr 2
    for gpcr_number, gpcr_sequence_number in zip(('GPCR1', 'GPCR2'), ('gpcr_sequence_number1', 'gpcr_sequence_number2')):
        # find which gpcr the position belongs to, to find the constructs
        if row[gpcr_number] == 'GIPR':
            constructs = gipr_constructs()
        elif row[gpcr_number] == 'GLP1R':
            constructs = glp1r_constructs()
        else:
            raise ValueError(f"Unknown GPCR: {row['GPCR1']}")
        gpcr = row[gpcr_number].lower()

        # for each construct, check if this sequence number matches the receptor
        for construct_number, construct_residues in constructs.items():
            if row[gpcr_sequence_number] in construct_residues:
                # if nan, replace
                constructs_here = df.at[i, 'Construct']
                if pd.isna(constructs_here):
                    df.at[i, 'Construct'] = construct_number
                    df.at[i, 'ConstructImpact'] = get_construct_impact()[gpcr][construct_number]
                else:
                    # if not nan, append
                    df.at[i, 'Construct'] = f"{constructs_here}|{construct_number}"
                    df.at[i, 'ConstructImpact'] = f"{df.at[i, 'ConstructImpact']}|{get_construct_impact()[gpcr][construct_number]}"

df.to_csv(out_dir / 'MI_inter_annotated_constructs.csv', index=False)