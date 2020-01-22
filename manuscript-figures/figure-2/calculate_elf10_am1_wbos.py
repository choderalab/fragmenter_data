"""
Filter small molecule of drugbank to include triple bonds and calculate WBOs. This is so that the mode at 3 shows up in
the distribution of WBOs

Ran on 2020-1-21
----------------
Software versions used:
fragmenter: 0.0.4+24.g18bdc15.dirty
openeye: 2019.Apr.2
"""

from fragmenter import chemi
from openeye import oechem, oedepict
import pandas as pd


def triple(smiles):
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, smiles)
    for b in mol.GetBonds():
        if b.IsInRing():
            continue
        if b.GetOrder() > 2:
            return True
    return False


def atomic_charge(smiles):
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, smiles)
    for a in mol.GetAtoms():
        if a.GetFormalCharge() != 0:
            return False
    return True


small_mols = pd.read_csv('drugbank_small_mols.csv')

small_mols['triple'] = small_mols['smiles'].apply(triple)
small_mols['charge'] = small_mols['smiles'].apply(atomic_charge)

filtered_drugbank = small_mols.loc[(small_mols['heavy_atoms'] <= 20) &
                                       (small_mols['fda_approved'] == True) &
                                       (small_mols['connected_components'] == 1) &
                                       (small_mols['charge'] == True) &
                                       (small_mols['triple'] == True),]

# Calculate WBOs
ofs = oechem.oemolostream('druglike_wbos.oeb')
ofs_2 = oechem.oemolostream('druglike_wbos.smi')
all_mols = []
for smiles in filtered_drugbank['smiles']:
    mol = chemi.smiles_to_oemol(smiles)
    all_mols.append(mol)
    try:
        charged_mol = chemi.get_charges(mol, strict_stereo=False, strict_types=False)
    except RuntimeError:
        print('charging failed for {}'.format(smiles))
    try:
        oechem.OEWriteMolecule(ofs, charged_mol)
        oechem.OEWriteMolecule(ofs_2, charged_mol)
    except Exception as e:
        print('Failed to charge {}'.format(smiles))
        print(e)


# Also include kinase inhibitors because they have a good amount of sulfur and phosphourous which have weaker bonds
mols = chemi.file_to_oemols('kinase_inhibitors.smi')
for mol in mols:
    try:
        charged_mol = chemi.get_charges(mol, strict_stereo=False, strict_types=False)
    except RuntimeError:
        print('charging failed for {}'.format(smiles))
    try:
        oechem.OEWriteMolecule(ofs, charged_mol)
        oechem.OEWriteMolecule(ofs_2, charged_mol)
    except Exception as e:
        print('Failed to charg')
all_mols.extend(mols)
