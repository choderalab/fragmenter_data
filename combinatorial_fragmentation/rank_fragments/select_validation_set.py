import json
import glob
from fragmenter import chemi
from openeye import oechem
import os
import shutil


# Collect all mmd scores
mols = glob.glob('validation_set/*/')
max_scores = []
failures = []
for name in mols:
    name = name.split('/')[1]
    print(name)
    try:
        with open('validation_set/{}/{}_mmd_exp_scores.json'.format(name, name), 'r') as f:
            scores = json.load(f)
    except FileNotFoundError:
        # Computation didn't complete
        failures.append(name)
        continue
    max_score = scores['greatest_discrepancy']
    for b in scores:
        if b == 'greatest_discrepancy':
            continue
        if max_score in scores[b]:
            bond = [int(i) for i in b.split('[')[-1].split(']')[0].split(',')]
    with open('validation_set/{}/{}_oe_wbo_with_score.json'.format(name, name), 'r') as f:
        frags = json.load(f)
        for bo in frags:
            for f in frags[bo]:
                if 'parent' in f:
                    parent_smiles = frags[bo][f]["map_to_parent"]
                   # parent_smiles = f.split('_')[0]
                    parent_mol = oechem.OEMol()
                    oechem.OESmilesToMol(parent_mol, parent_smiles)
                    # Get charge
                    charge = chemi.get_charge(parent_mol)
                    parent_mol.SetTitle(name)
                    break
            break

    max_scores.append((max_score, name, charge, parent_mol, bond))
print(len(max_scores))

# sort mol and visualize top 50 molecules.
sorted_scores = sorted(max_scores, reverse=True)
sorted_mols = [i[3] for i in sorted_scores[:100]]
sorted_bonds = [i[-1] for i in sorted_scores[:100]]
chemi.to_pdf(sorted_mols, fname='selected_validation_set.pdf',  bond_map_idx=sorted_bonds)

try:
    os.mkdir('selected')
except FileExistsError:
    print("overwriting selected")
validation_set = {}
for mol in sorted_scores[:100]:
    shutil.copytree('validation_set/{}'.format(mol[1]), 'selected/{}'.format(mol[1]))


