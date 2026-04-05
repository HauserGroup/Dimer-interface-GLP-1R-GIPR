reinitialize
bg_color white
set grid_mode, 1
load /Users/kcd635/HauserGroup Dropbox/luuk hd/GLP1R_GIPR_dimers/code/data/AF2_PDBs/glp1r_human___gipr_human.pdb, glp1_gipr
load /Users/kcd635/HauserGroup Dropbox/luuk hd/GLP1R_GIPR_dimers/code/data/AF2_PDBs/glp1r_human___gipr_human___glp1.pdb, glp1r_gipr_glp1
load /Users/kcd635/HauserGroup Dropbox/luuk hd/GLP1R_GIPR_dimers/code/data/AF2_PDBs/glp1r_human___gipr_human___gip.pdb, glp1r_gipr_gip

align glp1r_gipr_gip, glp1_gipr
align glp1r_gipr_glp1, glp1_gipr

# shared settings, selections and colors
@/Users/kcd635/HauserGroup Dropbox/luuk hd/GLP1R_GIPR_dimers/code/plots/pymol/commands/b_ray.pm
@/Users/kcd635/HauserGroup Dropbox/luuk hd/GLP1R_GIPR_dimers/code/plots/pymol/commands/a_selections.pm
@/Users/kcd635/HauserGroup Dropbox/luuk hd/GLP1R_GIPR_dimers/code/plots/pymol/commands/a_colors.pm

# show GPCRs as cartoon
show cartoon, glp1r
show cartoon, gipr
# show ligands as sticks
show sticks, peptide
hide cartoon, peptide
hide (not backbone)

# show interface residues as sticks
show sticks, glp1r_interface
show sticks, gipr_interface

# # flights
# show surface, glp1r_interface
# show surface, gipr_interface

# make sticks thicker, only backbone
set stick_radius, 0.3
# sticks for peptides, 0.5
set stick_radius, 0.5, peptide
hide (elem H)

# show surface of TMs
show surface, glp1r_tm4
show surface, gipr_tm1
show surface, gipr_tm2
set transparency, 0.6
set two_sided_lighting, on
set ray_shadows, off

color GLP1_blue, glp1r_tm4
color GIP_orange, gipr_tm1
color GIP_orange, gipr_tm2
color GLP1_blue, glp1r_interface
color GIP_orange, gipr_interface

# make cartoon less noticeable
cartoon loop
hide (elem H)
set cartoon_sampling, 10

# # label
# set label_size, 20
# label resi 58 and glp1r and name CA, 'GLP1R'
# label resi 54 and gipr and name CA, 'GIPR'
# label glp1r_gipr_glp1 and resi 28 and peptide and name HE2, 'GLP1'
# label glp1r_gipr_gip and resi 42 and peptide and name OE1, 'GIP'
set_view (        -0.803471148,   -0.154410169,    0.574965954,        -0.104786389,   -0.914020896,   -0.391896367,        0.586042702,   -0.375129342,    0.718210280,        0.000000000,    0.000000000, -447.346466064,        0.527187347,    3.873847961,    3.092056274,        -81704.851562500, 82599.523437500,  -20.000000000 )

#ray 3840, 2160
#png /Users/kcd635/Documents/GitHub/GLP1R_GIPR_dimers/plots/dimer_hypothesis_4k.png
