# tm4 goes from 261, 292
import pandas as pd

df = pd.read_csv(
    "/Users/kcd635/Documents/GitHub/GLP1R_GIPR_dimers/glp1r_human_chain_info.csv"
)
print(df)
df = df[df["protein_segment"] == "TM4"]

for i, row in df.iterrows():
    # aa = f'{row["amino_acid"]}{row["sequence_number"]}'
    aa = f'{row["sequence_number"]}'
    print(aa, end=",")
print()
