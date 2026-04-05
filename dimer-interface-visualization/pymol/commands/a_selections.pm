sele glp1r, chain A
deselect
sele gipr, chain B
deselect
sele peptide, chain C

# g275 from PNAS paper was not predicted to interact
sele glp1r_interface, resi 279+289+282+283+286+287+290+291+275 and glp1r
deselect
# R136 was included due to alpha-carbon proximity, but lacks side-chain interactions. 
sele gipr_interface, resi 136+147+192+150+151+146+142+143+188+139+191+138+141 and gipr
deselect

sele glp1r_tm4, resi 261-292 and glp1r
deselect
sele gipr_tm1, resi 124-162 and gipr
deselect
sele gipr_tm2, resi 167-194 and gipr
deselect

# constructs
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

# 1: 136 143 147 151
sele gipr_construct_1, resi 136+143+147+151 and gipr
deselect
# 2: 138 142 146 191
sele gipr_construct_2, resi 138+142+146+191 and gipr
deselect
# 3: 136 138 142 143 146 147 151 191
sele gipr_construct_3, resi 136+138+142+143+146+147+151+191 and gipr
deselect