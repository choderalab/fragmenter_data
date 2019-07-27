import cmiles
import openeye
import fragmenter

import json

import matplotlib.colors as mcolors



def add_atom_map(mol):
    """
    Add atom map to molecule while saving the 4-5 atom map (which specifies the central bond we are interested in
    drivine) in atom data.
    Parameters
    ----------
    mol : oemol
        Mol with 4, 5 atom map specifying bond of interest

    Returns
    -------

    """
    # copy molecule
    mol = openeye.oechem.OEMol(mol)
    if not cmiles.utils.has_explicit_hydrogen(mol):
        cmiles.utils.add_explicit_hydrogen(mol)
    cmiles._cmiles_oe.canonical_order_atoms(mol, in_place=True)
    for a in mol.GetAtoms():
        if a.GetMapIdx() != 0 :
            a.SetData('map_idx', a.GetMapIdx())
    cmiles.utils.remove_atom_map(mol, keep_map_data=True)
    # Add map indices on mol
    for a in mol.GetAtoms():
        a.SetMapIdx(a.GetIdx() + 1)
    return mol

def find_central_bond(mol):
    """
    Find the central bond of interest - it has 4 and 5 in the atom data
    Parameters
    ----------
    mol :

    Returns
    -------

    """
    central_bond = []
    for atom in mol.GetAtoms():
        if 'map_idx' in atom.GetData():
            if atom.GetData('map_idx') in [4, 5]:
                central_bond.append(atom.GetMapIdx())
    return tuple(central_bond)

fgroups =  [
    'dimethylamino',
    'methylamino',
    'amino',
    'ethylamino',
    'propylamino',
    'hydroxy',
    'methoxy',
    'ethoxy',
    'dimethylurea',
    'urea',
    'phenylurea',
    'ethylamide',
    'amide',
    'methyl',
    'carbamate',
    'benzoicacid',
    'ethoxycarbonyl',
    'trifluoromethyl',
    'trimethylamonium',
    'nitro']

color_keys = ['rosybrown', 'indianred', 'red', 'orange', 'gold', 'yellow','greenyellow', 'green', 'limegreen',
          'lightseagreen', 'teal', 'cyan', 'deepskyblue', 'royalblue', 'mediumslateblue', 'blueviolet', 'mediumorchid', 'violet',
          'palevioletred', 'lightpink']
colors = mcolors.CSS4_COLORS

# The indices here are slightly different than the indices already on QCArchive. This is to make sure we do not rerun them.
to_replace = {'amino': ('[H:4][NH:3][c:2]1[cH:1]cc(nc1)[N+](C)(C)C', '[H:4][NH:3][c:2]1ccc(n[cH:1]1)[N+](C)(C)C'),
            'amide': ('C[C:4](=O)[NH:3][c:2]1ccc(n[cH:1]1)[N+](C)(C)C', 'C[C:4](=O)[NH:3][c:2]1[cH:1]cc(nc1)[N+](C)(C)C'),
            'carbamate':
            ('C[N+](C)(C)c1c[c:2]([cH:1]nc1)[O:3][C:4](=O)N', 'C[N+](C)(C)c1[cH:1][c:2](cnc1)[O:3][C:4](=O)N'),
            'dimethylamino':
            ('C[N:3]([CH3:4])[c:2]1ccc(n[cH:1]1)[N+](C)(C)C', '[CH3:4][N:3](C)[c:2]1[cH:1]cc(nc1)[N+](C)(C)C'),
            'ethoxy':
            ('C[CH2:4][O:3][c:2]1ccc(n[cH:1]1)[N+](C)(C)C' ,'C[CH2:4][O:3][c:2]1[cH:1]cc(nc1)[N+](C)(C)C'),
            'ethylamino':
            ('C[CH2:4][NH:3][c:2]1ccc(n[cH:1]1)[N+](C)(C)C', 'C[CH2:4][NH:3][c:2]1[cH:1]cc(nc1)[N+](C)(C)C'),
            'methylamino':
            ('[CH3:4][NH:3][c:2]1ccc(n[cH:1]1)[N+](C)(C)C', '[CH3:4][NH:3][c:2]1[cH:1]cc(nc1)[N+](C)(C)C'),
            'phenylurea':
            ('C[N+](C)(C)c1c[cH:1][c:2](cn1)[NH:3][C:4](=O)N', 'C[N+](C)(C)c1c[cH:1][c:2](cn1)[NH:3][C:4](=O)N'),
            'urea':
            ('CN[C:4](=O)[NH:3][c:2]1ccc(n[cH:1]1)[N+](C)(C)C', 'CN[C:4](=O)[NH:3][c:2]1[cH:1]cc(nc1)[N+](C)(C)C')}

td_inputs = {}
for i, fgroup in enumerate(fgroups):
    print(fgroup)
    with open('data/{}_R1_wbos.json'.format(fgroup), 'r') as f:
        wbos = json.load(f)

    # Select molecules (up to 15) that are evenly distributed along WBOs
    wbos_dist = [round(wbo[0], 2) for wbo in wbos]
    differences = [j-i for i, j in zip(wbos_dist[:-1], wbos_dist[1:])]
    indices = [i for i, j in enumerate(differences) if j > 0]
    selected_mols = [wbos[i] for i in indices]
    selected_mols.append(wbos[-1])

    # Remove phenoxide because the negative charge is hard to deal with
    wbos_filtered = []
    for wbo, name, smiles in selected_mols:
        if not 'phenoxide' in name:
            wbos_filtered.append([wbo, name, smiles])
    if len(wbos_filtered) < 2:
        print('WBO range is too small. Dropping {}'.format(fgroup))
        continue
    if len(wbos_filtered) > 10:
        # remove some
        to_remove = wbos_filtered[1::int(len(wbos_filtered)/(len(wbos_filtered) - 10))]
        wbo_filtered_2 = []
        for j in wbos_filtered:
            if j in to_remove:
                continue
            else:
                wbo_filtered_2.append(j)
        wbo_filtered_2.append(wbos_filtered[-1])

        #wbo_filtered = list(filter(lambda a: a in to_remove, wbos_filtered))
        wbos_filtered = wbo_filtered_2

    smiles = [wbo[-1] for wbo in wbos_filtered]

    wbos = [wbo[0] for wbo in wbos_filtered]
    c_oe = openeye.oechem.OEColor(colors[color_keys[i]])
    fragmenter.chemi.to_pdf(smiles, fname='figures/qcarchive_torsiondrives/{}_to_drive.pdf'.format(fgroup), bo=wbos, bond_map_idx=(4, 5),
                 align=smiles[0], color=c_oe)


    print(fgroup, len(wbos_filtered))
    job_indices = []
    for wbo, name, smiles in wbos_filtered:
        mol =openeye.oechem.OEMol()
        openeye.oechem.OESmilesToMol(mol, smiles)
        mapped_mol = add_atom_map(mol)

        # Find bond
        central_bond = find_central_bond(mapped_mol)
        torsion = fragmenter.torsions.find_torsion_around_bond(mapped_mol, bond=central_bond)
        # get job index. Save job indices and WBOs for each fgroup - this will make analysis easier
        mapped_smiles = openeye.oechem.OEMolToSmiles(mapped_mol)
        job_index = cmiles.utils.to_canonical_label(openeye.oechem.OEMolToSmiles(mapped_mol), torsion)
        if fgroup in to_replace:

            if job_index == to_replace[fgroup][0]:
                print('replacing index for {}'.format(name))
                job_index = to_replace[fgroup][1]
        job_indices.append([job_index, torsion, name, wbo])

        # Generate conformers
        conformers = fragmenter.chemi.generate_conformers(mapped_mol, dense=True)
        cmiles_identifiers = cmiles.get_molecule_ids(mapped_smiles)

        # Check that mapping did not change
        if openeye.oechem.OEMolToSmiles(conformers) != mapped_smiles or mapped_smiles != cmiles_identifiers['canonical_isomeric_explicit_hydrogen_mapped_smiles']:
            print('Map changed for {}'.format(job_index))
        qcarchive_mols = [cmiles.utils.mol_to_map_ordered_qcschema(conformer, mapped_smiles) for conformer in conformers.GetConfs()]
        td_inputs[job_index] = {'input_molecules': qcarchive_mols,
                               'grid': [15],
                               'dihedral': torsion,
                               'cmiles_identifiers':  cmiles.get_molecule_ids(mapped_smiles),
                               'provenance': {'fragmenter_version': fragmenter.__version__,
                                              'openeye_version:': openeye.__version__,
                                              'dataset': {'fgroup': fgroup,
                                                          'name': name,
                                                          'oe_wbo': wbo}}
                               }

    with open('data/{}_td_job_indices.json'.format(fgroup), 'w') as f:
        json.dump(job_indices, f, indent=2, sort_keys=True)

with open('data/phenyl_set_torsiondrive_inputs.json', 'w') as f:
    json.dump(td_inputs, f, indent=2, sort_keys=True)

print('Total number of jobs')
print(len(td_inputs))

