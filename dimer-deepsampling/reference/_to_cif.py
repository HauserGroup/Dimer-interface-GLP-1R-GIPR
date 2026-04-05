from Bio.PDB import PDBParser, MMCIFIO
import pathlib
import os

if __name__ == "__main__":
    script_dir = pathlib.Path(__file__).parent.resolve()
    os.chdir(script_dir)
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("tmpl", "dimer_no_pep.pdb")
    io = MMCIFIO()
    io.set_structure(structure)
    io.save("dimer_no_pep.cif")
