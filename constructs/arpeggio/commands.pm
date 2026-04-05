reinitialize
load /Users/kcd635/Documents/GitHub/GLP1R_GIPR_dimers/constructs/arpeggio/xray-louisiana-charlie/xray-louisiana-charlie.pse

# set label size
set label_size, 20
set_name amide-amide, amide_amide
set_name amide-ring, amide_ring

# print(cmd.get_names('all'))
# ['upload_model_fixed', 'binding_site', 'undefined-proximal', 'undefined-vdw', 'polar-proximal', 'weakpolar-proximal', 'hbond-vdw', 'polar-vdw', 'ionic-vdwclash', 'polar-vdwclash', 'ionic-proximal', 'hbond-proximal', 'weakhbond-proximal', 'hbond-vdwclash', 'hydrophobic-proximal', 'aromatic-proximal', 'undefined-vdwclash', 'weakhbond-vdw', 'weakpolar-vdw', 'weakpolar-vdwclash', 'ionic-vdw', 'hydrophobic-vdw', 'aromatic-vdw', 'aromatic-vdwclash', 'hydrophobic-vdwclash', 'weakhbond-vdwclash', 'carbonyl-vdwclash', 'ring_centers', 'ring-EF', 'ring-OE', 'ring-OT', 'ring-FE', 'ring-ET', 'ring-EE', 'ring-OF', 'CARBONPI', 'METSULPHURPI', 'CATIONPI', 'DONORPI', 'amide-ring', 'amide_amide', 'sele']
hide dashes, (chain A not within 10 of chain B)
show dashes, (chain B not within 10 of chain A)
# dist interface_distances, (chain A), (chain B), 5

color wheat, chain B

# GLP1R chain A
# resi 279+287+283+286+279+291+290 = the interface of GLP1R 
hide sticks
show sticks, chain A and resi 279+287+283+286+279+291+290
# set opacity to 0.5 for chain A except resi 279+287+283+286+279+291+290
set cartoon_transparency, 0.8, chain A and not resi 279+287+283+286+279+291+290
# label every residue in the interface, only CA
label resi 279+287+283+286+279+291+290 and chain A and name CA, "%s %s" % (resn, resi)

# GIPR chain B
# resi 136+138+142+143+146+147+151+191 = the interface of GIPR
show sticks, chain B and resi 136+138+142+143+146+147+151+191
# set opacity to 0.5 for chain B except resi 136+138+142+143+146+147+151+191
set cartoon_transparency, 0.8, chain B and not resi 136+138+142+143+146+147+151+191
# label every residue in the interface, only CA
label resi 136+138+142+143+146+147+151+191 and chain B and name CA, "%s %s" % (resn, resi)

# add label connections to AA
set label_connector, 0
