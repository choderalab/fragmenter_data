import pandas as pd
from openeye import oechem
from openmoltools import openeye

# Import FDA approved kinase inhibitors (As of Jan 2018)
kinase_inhibitors = pd.read_csv('../clinical-kinase-inhibitors.csv')

# Generate mol from smiles input file
for i, inhibitor in kinase_inhibitors.iterrows():
    mol = openeye.smiles_to_oemol(inhibitor['smiles'])
    # generate conformers
    conformers = openeye.generate_conformers(mol)
    # write out mol2 files
    openeye.molecule_to_mol2(conformers, tripos_mol2_filename='{}.mol2'.format(inhibitor['inhibitor']), conformer=None)