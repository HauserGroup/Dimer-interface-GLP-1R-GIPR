from Bio import PDB
import argparse


def cif_to_pdb(cif_file, pdb_file):
    # Create a parser for CIF files
    parser = PDB.MMCIFParser(QUIET=True)

    # Read the structure from the CIF file
    structure = parser.get_structure("structure", cif_file)

    # Create a PDBIO object to save the structure as a PDB file
    io = PDB.PDBIO()
    io.set_structure(structure)

    # Write the structure to a PDB file
    io.save(pdb_file)


def run_main():
    # Create a command-line parser
    parser = argparse.ArgumentParser(description="Convert a CIF file to a PDB file")

    # Add the command-line arguments
    parser.add_argument("-i", type=str, help="Input CIF file")
    parser.add_argument("-o", type=str, help="Output PDB file")

    # Parse the command-line arguments
    args = parser.parse_args()

    # Convert the CIF file to a PDB file
    cif_to_pdb(args.i, args.o)


if __name__ == "__main__":
    run_main()
