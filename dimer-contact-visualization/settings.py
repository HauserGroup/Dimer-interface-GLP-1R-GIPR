def get_tm1_dict():
    "starts at 1, first amino acid is M1"
    return {"GIPR": (124, 162), "GLP1R": (138, 169)}


def get_protein_lens():
    "starts at 1"
    return {"GIPR": 466, "GLP1R": 463}


def get_tm4_dict():
    "starts at 1"
    return {"GIPR": (255, 282), "GLP1R": (261, 292)}


def get_chains():
    return {"GIPR": "B", "GLP1R": "A"}
