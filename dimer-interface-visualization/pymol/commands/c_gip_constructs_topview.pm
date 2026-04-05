@/Users/kcd635/HauserGroup Dropbox/luuk hd/GLP1R_GIPR_dimers/code/plots/pymol/commands/c_dimer_interface.pm

### cut below here and paste into script ###
set_view (\
    -0.217379034,    0.975525379,   -0.032991134,\
    -0.226198763,   -0.083226897,   -0.970515668,\
    -0.949512303,   -0.203509778,    0.238754287,\
    -5.000000000,    0.000000000,  -86.977516174,\
     3.869751453,   11.138881683,   -0.831684113,\
    -2.453570127,  176.408599854,  -20.000000000 )
### cut above here and paste into script ###

hide sticks
show sticks, gipr_construct_1
show sticks, gipr_construct_2

# Set the stick color to GIP_orange for the specific constructs
set stick_color, GIP_orange, gipr_construct_1
set stick_color, GIP_orange, gipr_construct_2

# Ensure the overall color for the cartoon remains as GIPR
color GIPR, gipr
color GIPR, gipr
color GLP1R, glp1r

### cut below here and paste into script ###
set_view (\
    -0.042215716,    0.981869578,    0.184747607,\
    -0.199293315,    0.172920197,   -0.964559793,\
    -0.979017556,   -0.077539757,    0.188385382,\
    -5.000000000,   -2.000000000,  -74.341949463,\
     3.869751453,   11.138881683,   -0.831684113,\
   -15.089080811,  163.773056030,  -20.000000000 )
### cut above here and paste into script ###

# test
set ray_trace_fog,0
set ray_shadows, 0
set depth_cue, 0
set cartoon_transparency, 0.4, GIPR
set cartoon_transparency, 0.4, GLP1R
####
show sticks, glp1r_interface
set stick_transparency, 0.2, glp1r_interface
hide (elem H)

##ray 3840, 2160
##png /Users/kcd635/HauserGroup Dropbox/luuk hd/GLP1R_GIPR_dimers/code/plots/dimer_gipr_constructs_topview_4k.png

# try again but don't make gipr trans
# set cartoon_transparency, 0.0
# set cartoon_transparency, 0.0, GIPR
# set cartoon_transparency, 0.0, GLP1R
# ray 3840, 2160
# png /Users/kcd635/HauserGroup Dropbox/luuk hd/GLP1R_GIPR_dimers/code/plots/dimer_gipr_constructs_topview_4k_non_transp.png