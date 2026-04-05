"""
1. match proteins by species.

"""
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pathlib
from Bio import SeqIO


def is_amino_acid_col(col, gpcr_name):
    if not col.startswith(gpcr_name):
        return False
    col = col.replace(f'{gpcr_name}_', '')
    if not col[0] in "ACDEFGHIKLMNPQRSTVWY":
        return False
    if not col[1:].isdigit():
        return False
    return True


def get_co_mutation_ratio(df, wt_rows, query_col, target_col):
    """ Calculate the co-mutation rate between two columns.

    from query_col -> target_col
    
    step 1:
        get WT residue from query_col
    """
    # print(f'wt aa is: {wt_rows}')
    # print(f'query col is {query_col}')
    wt_query_residue = wt_rows[query_col].values[0]
    wt_target_residue = wt_rows[target_col].values[0]
    
    # get species where the query is mutated
    wt_mutated_rows = df[df[query_col] != wt_query_residue]
    # print(f'Number of species with a mutated {query_col}: {len(wt_mutated_rows)}')
    if wt_mutated_rows.empty:
        return pd.NA
    
    # for this subset with a mutated query column, calculate the percent of co-mutated instances of the target cols
    co_mutated = wt_mutated_rows[wt_mutated_rows[target_col] != wt_target_residue]
    # print(f'Number of species with a co-mutated {target_col}: {len(co_mutated)}')
    co_mutated_ratio = len(co_mutated) / len(wt_mutated_rows)
    return co_mutated_ratio


if __name__ == "__main__":
    """
    """
    # get the species
    script_dir = pathlib.Path(__file__).parent
    df = pd.read_csv(script_dir / 'step4b_species_and_residues.csv')
    print(df)
    # species,id_glp1r,id_gipr,gipr_M1,gipr_T2,gipr_T3,gipr_S4,gipr_P5,gipr_I6,gipr_L7,gipr_Q8,gipr_L9,gipr_L10,gipr_L11,gipr_R12,gipr_L13,gipr_S14,gipr_L15,gipr_C16,gipr_G17,gipr_L18,gipr_L19,gipr_L20,gipr_Q21,gipr_R22,gipr_A23,gipr_E24,gipr_T25,gipr_G26,gipr_S27,gipr_K28,gipr_G29,gipr_Q30,gipr_T31,gipr_A32,gipr_G33,gipr_E34,gipr_L35,gipr_Y36,gipr_Q37,gipr_R38,gipr_W39,gipr_E40,gipr_R41,gipr_Y42,gipr_R43,gipr_R44,gipr_E45,gipr_C46,gipr_Q47,gipr_E48,gipr_T49,gipr_L50,gipr_A51,gipr_A52,gipr_A53,gipr_E54,gipr_P55,gipr_P56,gipr_S57,gipr_G58,gipr_L59,gipr_A60,gipr_C61,gipr_N62,gipr_G63,gipr_S64,gipr_F65,gipr_D66,gipr_M67,gipr_Y68,gipr_V69,gipr_C70,gipr_W71,gipr_D72,gipr_Y73,gipr_A74,gipr_A75,gipr_P76,gipr_N77,gipr_A78,gipr_T79,gipr_A80,gipr_R81,gipr_A82,gipr_S83,gipr_C84,gipr_P85,gipr_W86,gipr_Y87,gipr_L88,gipr_P89,gipr_W90,gipr_H91,gipr_H92,gipr_H93,gipr_V94,gipr_A95,gipr_A96,gipr_G97,gipr_F98,gipr_V99,gipr_L100,gipr_R101,gipr_Q102,gipr_C103,gipr_G104,gipr_S105,gipr_D106,gipr_G107,gipr_Q108,gipr_W109,gipr_G110,gipr_L111,gipr_W112,gipr_R113,gipr_D114,gipr_H115,gipr_T116,gipr_Q117,gipr_C118,gipr_E119,gipr_N120,gipr_P121,gipr_E122,gipr_K123,gipr_N124,gipr_E125,gipr_A126,gipr_F127,gipr_L128,gipr_D129,gipr_Q130,gipr_R131,gipr_L132,gipr_I133,gipr_L134,gipr_E135,gipr_R136,gipr_L137,gipr_Q138,gipr_V139,gipr_M140,gipr_Y141,gipr_T142,gipr_V143,gipr_G144,gipr_Y145,gipr_S146,gipr_L147,gipr_S148,gipr_L149,gipr_A150,gipr_T151,gipr_L152,gipr_L153,gipr_L154,gipr_A155,gipr_L156,gipr_L157,gipr_I158,gipr_L159,gipr_S160,gipr_L161,gipr_F162,gipr_R163,gipr_R164,gipr_L165,gipr_H166,gipr_C167,gipr_T168,gipr_R169,gipr_N170,gipr_Y171,gipr_I172,gipr_H173,gipr_I174,gipr_N175,gipr_L176,gipr_F177,gipr_T178,gipr_S179,gipr_F180,gipr_M181,gipr_L182,gipr_R183,gipr_A184,gipr_A185,gipr_A186,gipr_I187,gipr_L188,gipr_S189,gipr_R190,gipr_D191,gipr_R192,gipr_L193,gipr_L194,gipr_P195,gipr_R196,gipr_P197,gipr_G198,gipr_P199,gipr_Y200,gipr_L201,gipr_G202,gipr_D203,gipr_Q204,gipr_A205,gipr_L206,gipr_A207,gipr_L208,gipr_W209,gipr_N210,gipr_Q211,gipr_A212,gipr_L213,gipr_A214,gipr_A215,gipr_C216,gipr_R217,gipr_T218,gipr_A219,gipr_Q220,gipr_I221,gipr_V222,gipr_T223,gipr_Q224,gipr_Y225,gipr_C226,gipr_V227,gipr_G228,gipr_A229,gipr_N230,gipr_Y231,gipr_T232,gipr_W233,gipr_L234,gipr_L235,gipr_V236,gipr_E237,gipr_G238,gipr_V239,gipr_Y240,gipr_L241,gipr_H242,gipr_S243,gipr_L244,gipr_L245,gipr_V246,gipr_L247,gipr_V248,gipr_G249,gipr_G250,gipr_S251,gipr_E252,gipr_E253,gipr_G254,gipr_H255,gipr_F256,gipr_R257,gipr_Y258,gipr_Y259,gipr_L260,gipr_L261,gipr_L262,gipr_G263,gipr_W264,gipr_G265,gipr_A266,gipr_P267,gipr_A268,gipr_L269,gipr_F270,gipr_V271,gipr_I272,gipr_P273,gipr_W274,gipr_V275,gipr_I276,gipr_V277,gipr_R278,gipr_Y279,gipr_L280,gipr_Y281,gipr_E282,gipr_N283,gipr_T284,gipr_Q285,gipr_C286,gipr_W287,gipr_E288,gipr_R289,gipr_N290,gipr_E291,gipr_V292,gipr_K293,gipr_A294,gipr_I295,gipr_W296,gipr_W297,gipr_I298,gipr_I299,gipr_R300,gipr_T301,gipr_P302,gipr_I303,gipr_L304,gipr_M305,gipr_T306,gipr_I307,gipr_L308,gipr_I309,gipr_N310,gipr_F311,gipr_L312,gipr_I313,gipr_F314,gipr_I315,gipr_R316,gipr_I317,gipr_L318,gipr_G319,gipr_I320,gipr_L321,gipr_L322,gipr_S323,gipr_K324,gipr_L325,gipr_R326,gipr_T327,gipr_R328,gipr_Q329,gipr_M330,gipr_R331,gipr_C332,gipr_R333,gipr_D334,gipr_Y335,gipr_R336,gipr_L337,gipr_R338,gipr_L339,gipr_A340,gipr_R341,gipr_S342,gipr_T343,gipr_L344,gipr_T345,gipr_L346,gipr_V347,gipr_P348,gipr_L349,gipr_L350,gipr_G351,gipr_V352,gipr_H353,gipr_E354,gipr_V355,gipr_V356,gipr_F357,gipr_A358,gipr_P359,gipr_V360,gipr_T361,gipr_E362,gipr_E363,gipr_Q364,gipr_A365,gipr_R366,gipr_G367,gipr_A368,gipr_L369,gipr_R370,gipr_F371,gipr_A372,gipr_K373,gipr_L374,gipr_G375,gipr_F376,gipr_E377,gipr_I378,gipr_F379,gipr_L380,gipr_S381,gipr_S382,gipr_F383,gipr_Q384,gipr_G385,gipr_F386,gipr_L387,gipr_V388,gipr_S389,gipr_V390,gipr_L391,gipr_Y392,gipr_C393,gipr_F394,gipr_I395,gipr_N396,gipr_K397,gipr_E398,gipr_V399,gipr_Q400,gipr_S401,gipr_E402,gipr_I403,gipr_R404,gipr_R405,gipr_G406,gipr_W407,gipr_H408,gipr_H409,gipr_C410,gipr_R411,gipr_L412,gipr_R413,gipr_R414,gipr_S415,gipr_L416,gipr_G417,gipr_E418,gipr_E419,gipr_Q420,gipr_R421,gipr_Q422,gipr_L423,gipr_P424,gipr_E425,gipr_R426,gipr_A427,gipr_F428,gipr_R429,gipr_A430,gipr_L431,gipr_P432,gipr_S433,gipr_G434,gipr_S435,gipr_G436,gipr_P437,gipr_G438,gipr_E439,gipr_V440,gipr_P441,gipr_T442,gipr_S443,gipr_R444,gipr_G445,gipr_L446,gipr_S447,gipr_S448,gipr_G449,gipr_T450,gipr_L451,gipr_P452,gipr_G453,gipr_P454,gipr_G455,gipr_N456,gipr_E457,gipr_A458,gipr_S459,gipr_R460,gipr_E461,gipr_L462,gipr_E463,gipr_S464,gipr_Y465,gipr_C466,glp1r_M1,glp1r_A2,glp1r_G3,glp1r_A4,glp1r_P5,glp1r_G6,glp1r_P7,glp1r_L8,glp1r_R9,glp1r_L10,glp1r_A11,glp1r_L12,glp1r_L13,glp1r_L14,glp1r_L15,glp1r_G16,glp1r_M17,glp1r_V18,glp1r_G19,glp1r_R20,glp1r_A21,glp1r_G22,glp1r_P23,glp1r_R24,glp1r_P25,glp1r_Q26,glp1r_G27,glp1r_A28,glp1r_T29,glp1r_V30,glp1r_S31,glp1r_L32,glp1r_W33,glp1r_E34,glp1r_T35,glp1r_V36,glp1r_Q37,glp1r_K38,glp1r_W39,glp1r_R40,glp1r_E41,glp1r_Y42,glp1r_R43,glp1r_R44,glp1r_Q45,glp1r_C46,glp1r_Q47,glp1r_R48,glp1r_S49,glp1r_L50,glp1r_T51,glp1r_E52,glp1r_D53,glp1r_P54,glp1r_P55,glp1r_P56,glp1r_A57,glp1r_T58,glp1r_D59,glp1r_L60,glp1r_F61,glp1r_C62,glp1r_N63,glp1r_R64,glp1r_T65,glp1r_F66,glp1r_D67,glp1r_E68,glp1r_Y69,glp1r_A70,glp1r_C71,glp1r_W72,glp1r_P73,glp1r_D74,glp1r_G75,glp1r_E76,glp1r_P77,glp1r_G78,glp1r_S79,glp1r_F80,glp1r_V81,glp1r_N82,glp1r_V83,glp1r_S84,glp1r_C85,glp1r_P86,glp1r_W87,glp1r_Y88,glp1r_L89,glp1r_P90,glp1r_W91,glp1r_A92,glp1r_S93,glp1r_S94,glp1r_V95,glp1r_P96,glp1r_Q97,glp1r_G98,glp1r_H99,glp1r_V100,glp1r_Y101,glp1r_R102,glp1r_F103,glp1r_C104,glp1r_T105,glp1r_A106,glp1r_E107,glp1r_G108,glp1r_L109,glp1r_W110,glp1r_L111,glp1r_Q112,glp1r_K113,glp1r_D114,glp1r_N115,glp1r_S116,glp1r_S117,glp1r_L118,glp1r_P119,glp1r_W120,glp1r_R121,glp1r_D122,glp1r_L123,glp1r_S124,glp1r_E125,glp1r_C126,glp1r_E127,glp1r_E128,glp1r_S129,glp1r_K130,glp1r_R131,glp1r_G132,glp1r_E133,glp1r_R134,glp1r_S135,glp1r_S136,glp1r_P137,glp1r_E138,glp1r_E139,glp1r_Q140,glp1r_L141,glp1r_L142,glp1r_F143,glp1r_L144,glp1r_Y145,glp1r_I146,glp1r_I147,glp1r_Y148,glp1r_T149,glp1r_V150,glp1r_G151,glp1r_Y152,glp1r_A153,glp1r_L154,glp1r_S155,glp1r_F156,glp1r_S157,glp1r_A158,glp1r_L159,glp1r_V160,glp1r_I161,glp1r_A162,glp1r_S163,glp1r_A164,glp1r_I165,glp1r_L166,glp1r_L167,glp1r_G168,glp1r_F169,glp1r_R170,glp1r_H171,glp1r_L172,glp1r_H173,glp1r_C174,glp1r_T175,glp1r_R176,glp1r_N177,glp1r_Y178,glp1r_I179,glp1r_H180,glp1r_L181,glp1r_N182,glp1r_L183,glp1r_F184,glp1r_A185,glp1r_S186,glp1r_F187,glp1r_I188,glp1r_L189,glp1r_R190,glp1r_A191,glp1r_L192,glp1r_S193,glp1r_V194,glp1r_F195,glp1r_I196,glp1r_K197,glp1r_D198,glp1r_A199,glp1r_A200,glp1r_L201,glp1r_K202,glp1r_W203,glp1r_M204,glp1r_Y205,glp1r_S206,glp1r_T207,glp1r_A208,glp1r_A209,glp1r_Q210,glp1r_Q211,glp1r_H212,glp1r_Q213,glp1r_W214,glp1r_D215,glp1r_G216,glp1r_L217,glp1r_L218,glp1r_S219,glp1r_Y220,glp1r_Q221,glp1r_D222,glp1r_S223,glp1r_L224,glp1r_S225,glp1r_C226,glp1r_R227,glp1r_L228,glp1r_V229,glp1r_F230,glp1r_L231,glp1r_L232,glp1r_M233,glp1r_Q234,glp1r_Y235,glp1r_C236,glp1r_V237,glp1r_A238,glp1r_A239,glp1r_N240,glp1r_Y241,glp1r_Y242,glp1r_W243,glp1r_L244,glp1r_L245,glp1r_V246,glp1r_E247,glp1r_G248,glp1r_V249,glp1r_Y250,glp1r_L251,glp1r_Y252,glp1r_T253,glp1r_L254,glp1r_L255,glp1r_A256,glp1r_F257,glp1r_S258,glp1r_V259,glp1r_L260,glp1r_S261,glp1r_E262,glp1r_Q263,glp1r_W264,glp1r_I265,glp1r_F266,glp1r_R267,glp1r_L268,glp1r_Y269,glp1r_V270,glp1r_S271,glp1r_I272,glp1r_G273,glp1r_W274,glp1r_G275,glp1r_V276,glp1r_P277,glp1r_L278,glp1r_L279,glp1r_F280,glp1r_V281,glp1r_V282,glp1r_P283,glp1r_W284,glp1r_G285,glp1r_I286,glp1r_V287,glp1r_K288,glp1r_Y289,glp1r_L290,glp1r_Y291,glp1r_E292,glp1r_D293,glp1r_E294,glp1r_G295,glp1r_C296,glp1r_W297,glp1r_T298,glp1r_R299,glp1r_N300,glp1r_S301,glp1r_N302,glp1r_M303,glp1r_N304,glp1r_Y305,glp1r_W306,glp1r_L307,glp1r_I308,glp1r_I309,glp1r_R310,glp1r_L311,glp1r_P312,glp1r_I313,glp1r_L314,glp1r_F315,glp1r_A316,glp1r_I317,glp1r_G318,glp1r_V319,glp1r_N320,glp1r_F321,glp1r_L322,glp1r_I323,glp1r_F324,glp1r_V325,glp1r_R326,glp1r_V327,glp1r_I328,glp1r_C329,glp1r_I330,glp1r_V331,glp1r_V332,glp1r_S333,glp1r_K334,glp1r_L335,glp1r_K336,glp1r_A337,glp1r_N338,glp1r_L339,glp1r_M340,glp1r_C341,glp1r_K342,glp1r_T343,glp1r_D344,glp1r_I345,glp1r_K346,glp1r_C347,glp1r_R348,glp1r_L349,glp1r_A350,glp1r_K351,glp1r_S352,glp1r_T353,glp1r_L354,glp1r_T355,glp1r_L356,glp1r_I357,glp1r_P358,glp1r_L359,glp1r_L360,glp1r_G361,glp1r_T362,glp1r_H363,glp1r_E364,glp1r_V365,glp1r_I366,glp1r_F367,glp1r_A368,glp1r_F369,glp1r_V370,glp1r_M371,glp1r_D372,glp1r_E373,glp1r_H374,glp1r_A375,glp1r_R376,glp1r_G377,glp1r_T378,glp1r_L379,glp1r_R380,glp1r_F381,glp1r_I382,glp1r_K383,glp1r_L384,glp1r_F385,glp1r_T386,glp1r_E387,glp1r_L388,glp1r_S389,glp1r_F390,glp1r_T391,glp1r_S392,glp1r_F393,glp1r_Q394,glp1r_G395,glp1r_L396,glp1r_M397,glp1r_V398,glp1r_A399,glp1r_I400,glp1r_L401,glp1r_Y402,glp1r_C403,glp1r_F404,glp1r_V405,glp1r_N406,glp1r_N407,glp1r_E408,glp1r_V409,glp1r_Q410,glp1r_L411,glp1r_E412,glp1r_F413,glp1r_R414,glp1r_K415,glp1r_S416,glp1r_W417,glp1r_E418,glp1r_R419,glp1r_W420,glp1r_R421,glp1r_L422,glp1r_E423,glp1r_H424,glp1r_L425,glp1r_H426,glp1r_I427,glp1r_Q428,glp1r_R429,glp1r_D430,glp1r_S431,glp1r_S432,glp1r_M433,glp1r_K434,glp1r_P435,glp1r_L436,glp1r_K437,glp1r_C438,glp1r_P439,glp1r_T440,glp1r_S441,glp1r_S442,glp1r_L443,glp1r_S444,glp1r_S445,glp1r_G446,glp1r_A447,glp1r_T448,glp1r_A449,glp1r_G450,glp1r_S451,glp1r_S452,glp1r_M453,glp1r_Y454,glp1r_T455,glp1r_A456,glp1r_T457,glp1r_C458,glp1r_Q459,glp1r_A460,glp1r_S461,glp1r_C462,glp1r_S463
    # Castor canadensis,XP_020041914.1,XP_020028513.1,R,P,N,S,P,P,L,R,L,L,L,L,L,S,L,W,G,L,P,L,R,R,A,E,T,G,S,E,G,Q,T,T,R,A,L,Y,E,R,W,E,R,Y,R,Q,E,C,L,E,T,L,E,A,A,E,P,P,S,G,L,A,C,N,G,S,F,D,M,Y,V,C,W,D,S,T,A,P,N,A,T,A,R,A,S,C,P,W,Y,L,P,W,H,H,H,V,A,S,G,F,V,L,R,Q,C,G,S,D,G,Q,W,G,P,W,R,D,H,T,Q,C,E,N,P,E,K,N,G,A,F,Q,D,Q,R,L,I,L,E,W,L,Q,V,V,Y,T,I,G,Y,S,L,S,L,A,T,L,L,L,A,L,L,I,L,S,F,F,R,R,L,R,C,T,R,N,Y,I,H,M,N,L,F,T,S,F,M,L,R,A,A,A,I,L,T,R,D,R,L,L,S,A,P,G,P,Y,P,G,E,Q,T,R,T,P,W,N,Q,A,L,A,A,C,R,T,A,Q,I,V,T,Q,Y,C,V,G,A,N,Y,T,W,L,L,V,E,G,V,Y,L,H,R,L,L,V,L,V,G,G,S,E,D,G,H,F,H,C,Y,L,L,L,G,W,G,A,P,A,L,F,V,I,P,W,V,I,V,R,Y,L,Y,E,N,T,Q,C,W,E,R,N,D,V,K,A,I,W,W,I,I,R,T,P,I,L,I,T,I,L,I,N,F,L,I,F,M,R,I,L,G,I,L,V,S,K,L,R,T,R,Q,M,R,C,P,D,Y,R,L,R,L,A,R,S,T,L,T,L,V,P,L,L,G,V,H,E,V,V,F,A,P,V,T,E,E,Q,A,Q,G,A,L,R,F,T,K,L,G,F,E,I,F,L,R,S,F,Q,G,F,L,V,S,V,L,Y,C,F,V,N,K,E,V,Q,S,E,I,R,R,G,W,H,H,R,R,L,R,R,S,L,S,E,E,S,R,Q,P,P,E,H,A,S,R,A,G,P,P,S,C,V,P,G,E,V,P,I,G,C,A,L,T,S,G,A,L,P,G,R,G,G,E,G,S,R,V,L,E,S,Y,C,M,A,G,A,P,S,A,L,R,L,A,L,L,L,L,G,A,V,G,R,A,G,P,R,P,Q,G,A,T,V,S,L,S,E,T,V,Q,K,W,R,E,Y,R,R,Q,C,Q,R,F,L,I,E,A,P,P,P,A,T,G,L,F,C,N,R,T,F,D,D,Y,A,C,W,P,D,G,S,P,G,S,F,V,N,V,S,C,P,W,Y,L,P,W,A,S,S,V,L,Q,G,H,V,Y,R,F,C,T,A,E,G,L,W,L,H,K,D,N,S,S,L,P,W,R,D,L,S,E,C,E,E,S,K,R,G,E,K,S,S,L,E,E,Q,L,L,S,L,Y,V,I,Y,T,V,G,Y,A,L,S,F,S,A,L,V,I,A,S,A,I,L,L,S,F,R,H,L,H,C,T,R,N,Y,I,H,L,N,L,F,A,S,F,I,L,R,A,L,S,V,F,I,K,D,A,V,L,K,W,M,Y,S,T,A,A,Q,Q,H,Q,W,D,G,L,L,S,Y,Q,D,S,L,G,C,R,L,V,F,L,L,M,Q,Y,C,V,A,A,N,Y,Y,W,L,L,V,E,G,V,Y,L,Y,T,L,L,A,F,S,V,F,S,E,Q,R,I,F,R,L,Y,L,S,I,G,W,G,V,P,L,L,F,V,I,P,W,G,I,V,K,Y,L,Y,E,D,E,G,C,W,T,R,N,S,N,M,N,Y,W,L,I,I,R,L,P,I,L,F,A,I,G,V,N,F,L,I,F,I,R,V,I,C,I,V,V,S,K,L,K,A,N,L,M,C,K,T,D,I,K,C,R,L,A,K,S,T,L,T,L,I,P,L,L,G,T,H,E,V,I,F,A,F,V,M,D,E,H,A,R,G,T,L,R,F,I,K,L,F,T,E,L,S,F,T,S,F,Q,G,L,M,V,A,I,L,Y,C,F,V,N,N,E,V,Q,M,E,F,R,K,S,W,E,R,W,R,L,E,H,L,H,I,Q,R,D,S,S,M,K,P,L,K,C,P,T,S,S,L,S,S,G,G,T,V,G,S,S,V,Y,A,A,T,C,Q,A,S,C,S

    # per amino acid (gipr_M1 ... glp1r_S463), get the co-mutation per residue
    glp1r_cols = [c for c in df.columns if is_amino_acid_col(c, 'glp1r')]
    gipr_cols = [c for c in df.columns if is_amino_acid_col(c, 'gipr')]
   
    # get the co-mutation empty matrix
    glp1r_co_mutations = pd.DataFrame(index=glp1r_cols, columns=gipr_cols, data=pd.NA)
    
    # fill using get_co_mutation_rate
    for glp1r_col in glp1r_cols:
        print(f'Processing {glp1r_col}')
        wt_rows = df[df['species'] == 'Homo sapiens']
        # apply get_co_mutation_ratio to all gipr_cols
        for gipr_col in gipr_cols:
            glp1r_co_mutations.loc[glp1r_col, gipr_col] = get_co_mutation_ratio(df, wt_rows, glp1r_col, gipr_col)

    glp1r_co_mutations.to_csv(script_dir / 'step5_glp1r_gipr_co_mutations.csv')
