"""
Visualize for lifemap
"""
import pandas as pd 
import pathlib

script_dir = pathlib.Path(__file__).parent.resolve()
df = pd.read_csv(script_dir / 'step4b_species_and_residues.csv')
species_names = df['species'].unique()
with open(script_dir / 'step4c_species_names.txt', 'w') as f:
    for name in species_names:
        f.write(f"{name}\n")

input("Species names written to step4c_species_names.txt. Upload to https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi and create step4c_species_names_tax_report.txt...")

# Read the tax report
tax_report = pd.read_csv(script_dir / 'step4c_species_names_tax_report.txt', sep='\t')
tax_report.columns = tax_report.columns.str.strip()
taxids = tax_report['taxid'].unique()
# write taxids
with open(script_dir / 'step4c_taxids.txt', 'w') as f:
    for taxid in taxids:
        f.write(f"{taxid} ")
