reinitialize
bg_color white

# load
load /Users/kcd635/HauserGroup Dropbox/luuk hd/GLP1R_GIPR_dimers/code/data/AF2_PDBs/glp1r_human___gipr_human.pdb, glp1_gipr
load /Users/kcd635/HauserGroup Dropbox/luuk hd/GLP1R_GIPR_dimers/code/data/AF2_PDBs/glp1r_human___gipr_human___glp1.pdb, glp1r_gipr_glp1
load /Users/kcd635/HauserGroup Dropbox/luuk hd/GLP1R_GIPR_dimers/code/data/AF2_PDBs/glp1r_human___gipr_human___gip.pdb, glp1r_gipr_gip
align glp1r_gipr_gip, glp1_gipr
align glp1r_gipr_glp1, glp1_gipr

# delete all except the dimer
delete not glp1_gipr

# shared settings, selections and colors
@/Users/kcd635/HauserGroup Dropbox/luuk hd/GLP1R_GIPR_dimers/code/plots/pymol/commands/b_ray.pm
@/Users/kcd635/HauserGroup Dropbox/luuk hd/GLP1R_GIPR_dimers/code/plots/pymol/commands/a_selections.pm
@/Users/kcd635/HauserGroup Dropbox/luuk hd/GLP1R_GIPR_dimers/code/plots/pymol/commands/a_colors.pm

hide cartoon, peptide

# make the receptors white and transparent cartoons
hide cartoon, glp1r 
hide cartoon, gipr
# color white, glp1r
# color white, gipr

show cartoon, glp1r_tm4
show cartoon, gipr_tm1
show cartoon, gipr_tm2

set cartoon_transparency, 0.5, not glp1r_interface and not gipr_interface

color GLP1_blue, glp1r_interface
color GIP_orange, gipr_interface
show sticks, glp1r_interface
show sticks, gipr_interface

# make sticks thicker, only backbone
set stick_radius, 0.2

# make cartoon less noticeable
cartoon loop
hide (elem H)

# set cartoon_oval_length , 0.2


set_view (        -0.803471148,   -0.154410169,    0.574965954,        -0.104786389,   -0.914020896,   -0.391896367,        0.586042702,   -0.375129342,    0.718210280,        0.000000000,    0.000000000, -447.346466064,        0.527187347,    3.873847961,    3.092056274,        -81704.851562500, 82599.523437500,  -20.000000000 )

### cut below here and paste into script ###
set_view (\
    -0.803471148,   -0.154410169,    0.574965954,\
    -0.104786389,   -0.914020896,   -0.391896367,\
     0.586042702,   -0.375129342,    0.718210280,\
     0.000000000,   -0.000000000, -197.172775269,\
     3.869751453,   11.138881683,   -0.831684113,\
   155.452514648,  238.893035889,  -20.000000000 )
### cut above here and paste into script ###
###ray 3840, 2160; png /Users/kcd635/HauserGroup Dropbox/luuk hd/GLP1R_GIPR_dimers/code/plots/dimer_interface_4k.png


### cut below here and paste into script ###
set_view (\
    -0.787376940,    0.300336152,   -0.538361192,\
     0.111520231,   -0.789496303,   -0.603540421,\
    -0.606296778,   -0.535252333,    0.588135183,\
     0.000000000,    0.000000000, -130.699386597,\
     3.869751453,   11.138881683,   -0.831684113,\
    88.979148865,  172.419631958,  -20.000000000 )
### cut above here and paste into script ###
###ray 3840, 2160
###png /Users/kcd635/HauserGroup Dropbox/luuk hd/GLP1R_GIPR_dimers/code/plots/dimer_interface_4k_rotated.png


# add labels for each interface residue
label glp1r_interface and name CA, "%s" % resi
label gipr_interface and name CA, "%s" % resi
# labels -1 = 1 angst
set label_size, -1.25

# set label_color, yellow 
# set label_bg_color, black
set label_outline_color, black
set label_position, (0, 0, 1)
