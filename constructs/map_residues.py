import pathlib


residue_map = {
    "V259A": "V282A",
    "L256A": "L279A",
    "G252A": "G275A",  # GvplLfv
}

# mapped the positives from the paper to residues from uniprot
positives = [
    # ["L256A", "V259A"],
    ["L279A", "V282A"],
    # ["G252A", "L256A", "V259A"],
    ["G275A", "L279A", "V282A"],
]

# the negatives (dimer or trimer mutations on TM1 (turned away from the interface))
negatives = [
    ["F143A", "Y145A"],
    ["F143A", "I146A"],
    ["Y145A", "Y148A"],
    ["Y148A", "G151A"],
    ["T149A", "Y152A"],
    # =======
    ["E139A", "F143A", "Y145A"],
    ["Q140A", "F143A", "I146A"],
    ["F143A", "Y145A", "Y148A"],
    ["Y145A", "Y148A", "G151A"],
    ["I146A", "T149A", "Y152A"],
]


glp1r_seq = "MAGAPGPLRLALLLLGMVGRAGPRPQGATVSLWETVQKWREYRRQCQRSLTEDPPPATDLFCNRTFDEYACWPDGEPGSFVNVSCPWYLPWASSVPQGHVYRFCTAEGLWLQKDNSSLPWRDLSECEESKRGERSSPEEQLLFLYIIYTVGYALSFSALVIASAILLGFRHLHCTRNYIHLNLFASFILRALSVFIKDAALKWMYSTAAQQHQWDGLLSYQDSLSCRLVFLLMQYCVAANYYWLLVEGVYLYTLLAFSVLSEQWIFRLYVSIGWGVPLLFVVPWGIVKYLYEDEGCWTRNSNMNYWLIIRLPILFAIGVNFLIFVRVICIVVSKLKANLMCKTDIKCRLAKSTLTLIPLLGTHEVIFAFVMDEHARGTLRFIKLFTELSFTSFQGLMVAILYCFVNNEVQLEFRKSWERWRLEHLHIQRDSSMKPLKCPTSSLSSGATAGSSMYTATCQASCS"


def create_mutated_fasta(
    wt_seq: str, mutations: list, save_p: pathlib.Path, header=">seq"
) -> str:
    """Create a fasta file that applies the mutations to the wildtype sequence."""
    wt_seq = wt_seq.upper()
    mut_seq = wt_seq
    for mutation in mutations:
        mutation = mutation.upper()
        wt_res = mutation[0]
        mut_res = mutation[-1]
        mut_index = int(mutation[1:-1]) - 1

        # assert that the wildtype residue is correct
        assert wt_seq[mut_index] == wt_res

        # apply the mutation
        mut_seq = mut_seq[:mut_index] + mut_res + mut_seq[mut_index + 1 :]

    # s ave the fasta file
    with open(save_p, "w") as f:
        f.write(f"{header}\n{mut_seq}\n")


def main():
    script_dir = pathlib.Path(__file__).resolve().parent
    script_dir = script_dir / "InverseFold"
    fasta_dir = script_dir / "fastas"
    fasta_dir.mkdir(exist_ok=True)

    # create mutations for the positives
    pos_dir = fasta_dir / "positives"
    pos_dir.mkdir(exist_ok=True)
    for mutation_list in positives:
        header = "_".join(mutation_list)
        save_p = pos_dir / f"{header}.fasta"
        header = f">{header}"
        create_mutated_fasta(glp1r_seq, mutation_list, save_p, header=header)

    # create mutations for the negatives
    neg_dir = fasta_dir / "negatives"
    neg_dir.mkdir(exist_ok=True)
    for mutation_list in negatives:
        header = "_".join(mutation_list)
        save_p = neg_dir / f"{header}.fasta"
        header = f">{header}"
        create_mutated_fasta(glp1r_seq, mutation_list, save_p, header=header)


def create_mutation_lists_mcsm():
    """
    Provide one mutation per line

    Format:
    1st column: Chain ID
    2nd column: Mutation (wild-type, residue position and mutant)

    Example:
    A R282W
    A Y236F
    A N239Y
    A C242S
    """
    script_dir = pathlib.Path(__file__).resolve().parent
    out_dir = script_dir / "mcsm/positives"
    out_dir.mkdir(exist_ok=True)
    for mutation_list in positives:
        header = "_".join(mutation_list)
        save_p = out_dir / f"{header}.txt"
        with open(save_p, "w") as f:
            for mutation in mutation_list:
                f.write(f"A {mutation}\n")

    out_dir = script_dir / "mcsm/negatives"
    out_dir.mkdir(exist_ok=True)
    for mutation_list in negatives:
        header = "_".join(mutation_list)
        save_p = out_dir / f"{header}.txt"
        with open(save_p, "w") as f:
            for mutation in mutation_list:
                f.write(f"A {mutation}\n")


if __name__ == "__main__":
    # ESMFold
    main()

    # MCSM
    create_mutation_lists_mcsm()
