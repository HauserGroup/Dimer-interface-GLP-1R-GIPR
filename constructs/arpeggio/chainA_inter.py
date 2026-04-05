"""
"""
import pathlib


def main():
    script_dir = pathlib.Path(__file__).parent.absolute()
    residues = script_dir / 'chainA_inter_5ang.txt'
    with open(residues, 'r') as f:
        lines = f.readlines()

    for line in lines:
        amino_acid_pos, amino_acid = line.split()
        print(f'/A/{amino_acid_pos}/')

if __name__ == '__main__':
    main()
