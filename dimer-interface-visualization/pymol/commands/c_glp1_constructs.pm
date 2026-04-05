@/Users/kcd635/HauserGroup Dropbox/luuk hd/GLP1R_GIPR_dimers/code/plots/pymol/commands/c_dimer_interface.pm


### cut below here and paste into script ###
set_view (\
    -0.981306493,    0.136465207,    0.135665610,\
    -0.192425415,   -0.698825359,   -0.688919246,\
     0.000793331,   -0.702150106,    0.712019622,\
     0.000000000,   -2.000000000,  -90.800354004,\
     3.869751453,   11.138881683,   -0.831684113,\
    49.080135345,  132.520584106,  -20.000000000 )
### cut above here and paste into script ###


deselect

# Ensure the overall color for the cartoon remains as GIPR
color GIPR, gipr
color GLP1R, glp1r


# hide cartoon
show cartoon, glp1r_tm4 and gipr_tm1 and gipr_tm2
show sticks, glp1r_construct_1 and glp1r_construct_2 and glp1r_construct_3 and glp1r_construct_4
hide sticks, gipr_interface
hide sticks, glp1r_interface

# test
set ray_trace_fog,0
set ray_shadows, 0
set depth_cue, 0
set cartoon_transparency, 0.5, GIPR
set cartoon_transparency, 0.5, GLP1R

# construct 1
set stick_color, GLP1R_blue, glp1r_construct_2
set stick_transparency, 0.5, glp1r_construct_2
set stick_color, GLP1R_blue, glp1r_construct_3
set stick_transparency, 0.5, glp1r_construct_3
set stick_color, GLP1_blue, glp1r_construct_1
hide sticks, glp1r_construct_2
hide sticks, glp1r_construct_3
show sticks, glp1r_construct_1
set stick_transparency, 0.0, glp1r_construct_1
hide (elem H)
##ray 3840, 2160
##png /Users/kcd635/HauserGroup Dropbox/luuk hd/GLP1R_GIPR_dimers/code/plots/dimer_glp1r_construct1_4k.png

# # construct 2
set stick_color, GLP1R_blue, glp1r_construct_1
set stick_transparency, 0.5, glp1r_construct_1
set stick_color, GLP1R_blue, glp1r_construct_3
set stick_transparency, 0.5, glp1r_construct_3
set stick_color, GLP1_blue, glp1r_construct_2
hide sticks, glp1r_construct_1
hide sticks, glp1r_construct_3
show sticks, glp1r_construct_2
set stick_transparency, 0.0, glp1r_construct_2
hide (elem H)
##ray 3840, 2160
##png /Users/kcd635/HauserGroup Dropbox/luuk hd/GLP1R_GIPR_dimers/code/plots/dimer_glp1r_construct2_4k.png

# # construct 3
set stick_color, GLP1R_blue, glp1r_construct_1
set stick_transparency, 0.5, glp1r_construct_1
set stick_color, GLP1R_blue, glp1r_construct_2
set stick_transparency, 0.5, glp1r_construct_2
set stick_color, GLP1_blue, glp1r_construct_3
set stick_transparency, 0.0, glp1r_construct_3
hide sticks, glp1r_construct_1
hide sticks, glp1r_construct_2
show sticks, glp1r_construct_3
hide (elem H)
##ray 3840, 2160
##png /Users/kcd635/HauserGroup Dropbox/luuk hd/GLP1R_GIPR_dimers/code/plots/dimer_glp1r_construct3_4k.png

