# COLOUR BY RESIDUE - PATRICK BRENNAN, UNIVERSITY OF OXFORD - AUGUST 2020

from pymol import cmd


def ngl(selection="all"):
    cmd.select("backbone_carbons", "name C+CA")
    cmd.select("ALA", "name CB and resn ALA")
    cmd.select("ARG", "name CB+CG+CD+CZ and resn ARG")
    cmd.select("ASN", "name CB+CG and resn ASN")
    cmd.select("ASP", "name CB+CG and resn ASP")
    cmd.select("CYS", "name CB and resn CYS")
    cmd.select("GLN", "name CB+CG+CD and resn GLN")
    cmd.select("GLU", "name CB+CG+CD and resn GLU")
    cmd.select("HIS", "name CB+CG+CD2+CE1 and resn HIS")
    cmd.select("ILE", "name CB+CG1+CG2+CD1 and resn ILE")
    cmd.select("LEU", "name CB+CG+CD1+CD2 and resn LEU")
    cmd.select("LYS", "name CB+CG+CD+CE and resn LYS")
    cmd.select("MET", "name CB+CG+CE and resn MET")
    cmd.select("PHE", "name CB+CG+CD1+CD2+CE1+CE2+CZ and resn PHE")
    cmd.select("PRO", "name CB+CG+CD and resn PRO")
    cmd.select("SER", "name CB and resn SER")
    cmd.select("THR", "name CB+CG2 and resn THR")
    cmd.select("TRP", "name CB+CG+CD1+CD2+CE2+CE3+CZ2+CZ3+CH2 and resn TRP")
    cmd.select("TYR", "name CB+CG+CD1+CD2+CE1+CE2+CZ and resn TYR")
    cmd.select("VAL", "name CB+CG1+CG2 and resn VAL")

    code = {
        "backbone_carbons": "white",
        "ALA": "lime",
        "ARG": "density",
        "ASN": "deepsalmon",
        "ASP": "warmpink",
        "CYS": "paleyellow",
        "GLN": "tv_red",
        "GLU": "ruby",
        "HIS": "slate",
        "ILE": "forest",
        "LEU": "smudge",
        "LYS": "deepblue",
        "MET": "sand",
        "PHE": "gray40",
        "PRO": "gray20",
        "SER": "tv_orange",
        "THR": "brown",
        "TRP": "palegreen",
        "TYR": "wheat",
        "VAL": "pink",
    }
    cmd.select("none")
    for elem in code:
        line = "color " + code[elem] + "," + elem + "&" + selection
        cmd.do(line)


cmd.extend("ngl", ngl)
