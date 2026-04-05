from pymol import cmd

def print_contacts(residue_selection, other_chain_selection, distance_cutoff=4.0):
    # Identify contacts within the distance cutoff between the residue and the other chain
    contacts = cmd.find_pairs(residue_selection, other_chain_selection, cutoff=distance_cutoff)

    # Print the contacts and interaction types
    for (atom1, atom2) in contacts:
        # Get coordinates for each atom in the contact
        res1_info = cmd.get_model(atom1).atom[0]
        res2_info = cmd.get_model(atom2).atom[0]
        
        # Calculate the distance manually
        coords1 = cmd.get_atom_coords(atom1)
        coords2 = cmd.get_atom_coords(atom2)
        distance = ((coords1[0] - coords2[0])**2 + (coords1[1] - coords2[1])**2 + (coords1[2] - coords2[2])**2) ** 0.5

        # Define interaction type based on distance (customize thresholds as needed)
        if distance < 3.5:
            interaction_type = "Hydrogen bond" if "O" in (res1_info.name, res2_info.name) or "N" in (res1_info.name, res2_info.name) else "Van der Waals"
        elif distance < 4.0:
            interaction_type = "Van der Waals"
        else:
            interaction_type = "Other"

        # Print information
        print(f"Residue {res1_info.resn}{res1_info.resi} ({res1_info.chain}) - "
              f"Residue {res2_info.resn}{res2_info.resi} ({res2_info.chain}) | "
              f"Distance: {distance:.2f} Å | Interaction Type: {interaction_type}")


# Example usage:
# Select residue of interest and other chain, e.g., "resi 100 and chain A" and "chain B"
DIMER = 'upload_model_fixed'
GLP1R_CHAIN = "A"
GIPR_CHAIN = "B"


# residues of interest
# GIPR Q138, T142, T146, D191 
for RES in ["138", "142", "146", "191"]:
    print_contacts(f"{DIMER} and resi {RES} and chain {GIPR_CHAIN}", f"{DIMER} and chain {GLP1R_CHAIN}", distance_cutoff=4.0)

print("\n\n\n\n")
# residues of interest
# GLP1R G275, L279, and V282
for RES in ["275", "279", "282"]:
    print_contacts(f"{DIMER} and resi {RES} and chain {GLP1R_CHAIN}", f"{DIMER} and chain {GIPR_CHAIN}", distance_cutoff=4.0)
    