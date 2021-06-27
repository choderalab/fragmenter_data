"""
OFF default level of theory (B3LYP-D3(BJ) / DZVP) is expensive for size of kinase inhibitor molecules.
Since we are trying to see how good/bad hf3c is relative to AM1, we don't need to calculate WBO distributions
for all kinase ihibitors. This script filters out inhibitors that have less than 30 heavy atoms. It also adds
Imatinib and Gefitinib even if they are bigger because they are examples used in the paper
"""

import json
import gzip
from openeye import oechem

with gzip.open('optimization_inputs.json.gz', 'r') as f:
    mols = json.load(f)

smiles = []
heavy_atoms = []
for m in mols:
    smile = m['cmiles_identifiers']['canonical_isomeric_smiles']
    oemol = oechem.OEMol()
    oechem.OESmilesToMol(oemol, smile)
    n = 0
    for a in oemol.GetAtoms():
        if not a.IsHydrogen():
            n += 1
    smiles.append(smile)
    heavy_atoms.append(n)

to_compute = ['Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(CC5)C',
              'COc1cc2c(cc1OCCCN3CCOCC3)c(ncn2)Nc4ccc(c(c4)Cl)F']

indices = []
for i, n in enumerate(heavy_atoms):
    if smiles[i] in to_compute:
        indices.append(smiles[i])
    if n < 30:
        indices.append(smiles[i])

all_indices = []
for mol in mols:
    smile = mol['cmiles_identifiers']['canonical_isomeric_smiles']
    if smile in indices:
        for i, conf in enumerate(mol['initial_molecules']):
            all_indices.append('{}-{}'.format(smile, str(i)))

print(len(all_indices))

with open('indices_to_comput.json', 'w') as f:
    json.dump(all_indices, f, indent=2, sort_keys=True)