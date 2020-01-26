import fragmenter
import cmiles
from openeye import oechem
import json

# Get cmiles identifiers for kinase inhibitors
mols = fragmenter.chemi.file_to_oemols('kinase_inhibitors.smi')
cmiles_identifiers = {}
for mol in mols:
    name = mol.GetTitle()
    cmiles_identifiers[name] = cmiles.get_molecule_ids(oechem.OEMolToSmiles(mol), strict=False)
with open('data/kinase_inhibitors_cmiles_ids.json', 'w') as f:
    json.dump(cmiles_identifiers, f, indent=2, sort_keys=True)