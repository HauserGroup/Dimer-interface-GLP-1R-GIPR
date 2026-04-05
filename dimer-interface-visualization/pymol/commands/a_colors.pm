set_color GLP1R_blue = [0.843, 0.922, 1.000]
set_color GIPR_orange = [1.0, 0.871, 0.706]
set_color GLP1_blue = [0.329, 0.627, 0.988]
set_color GIP_orange = [1.0, 0.486, 0.008]



# always color chains the same
# 215, 235, 255
color GLP1R_blue, chain A
# 255, 222, 180
color GIPR_orange, chain B
# 84, 160, 252
color GLP1_blue, chain C and glp1r_gipr_glp1
# 255, 124, 2
color GIP_orange, chain C and glp1r_gipr_gip

# gip 1 construct is lighter. 255 224 193
set_color gip_construct_1_color = [1.0, 0.878, 0.757]
# gip 2 construct is darker. 245 195 137
set_color gip_construct_2_color = [0.961, 0.765, 0.537]