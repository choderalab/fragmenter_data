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

# Do not drive trimethylamonium and trifluoromethyl because it is too congested
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
    'nitro']

color_keys = ['rosybrown', 'indianred', 'red', 'orange', 'gold', 'yellow','greenyellow', 'green', 'limegreen',
          'lightseagreen', 'teal', 'cyan', 'deepskyblue', 'royalblue', 'mediumslateblue', 'blueviolet', 'mediumorchid',
           'lightpink']
colors = mcolors.CSS4_COLORS

# The indices here are slightly different than the indices already on QCArchive. This is to make sure we do not rerun them.
# Added later: In the end we did not need to replace because this set was run with the default OFF specs so all jobs
# were rerun. But the job indices are still as in this dictionary so it is here for the record.
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
            [('C[N+](C)(C)c1cc[c:2]([cH:1]n1)[NH:3][C:4](=O)N', 'C[N+](C)(C)c1c[cH:1][c:2](cn1)[NH:3][C:4](=O)N'),
            ('c1cc([cH:1][c:2](c1)[NH:3][C:4](=O)N)C(=O)O', 'c1cc(c[c:2]([cH:1]1)[NH:3][C:4](=O)N)C(=O)O')],
            'urea':
            ('CN[C:4](=O)[NH:3][c:2]1ccc(n[cH:1]1)[N+](C)(C)C', 'CN[C:4](=O)[NH:3][c:2]1[cH:1]cc(nc1)[N+](C)(C)C'),
            'ethoxycarbonyl':
            ('CC[O:4][C:3](=O)[c:2]1ccnc([cH:1]1)[N+](C)(C)C', 'CCO[C:3](=[O:4])[c:2]1[cH:1]cnc(c1)[N+](C)(C)C'),
            'ethylamide':
            ( '[H:4][N:3]([c:2]1ccc(c[cH:1]1)F)C(=O)CC','CC[C:4](=O)[NH:3][c:2]1[cH:1]cc(cc1)F'),
             'nitro':
             ('C[N+](C)(C)c1c[c:2]([cH:1]cn1)[N+:3](=[O:4])[O-]', 'C[N+](C)(C)c1c[c:2]([cH:1]cn1)[N+:3](=O)[O-:4]'),
             'benzoicacid':
             ('C[N+](C)(C)c1[cH:1][c:2](ccn1)[C:3](=O)[OH:4]', 'C[N+](C)(C)c1c[c:2]([cH:1]cn1)[C:3](=[O:4])O')
            }
td_inputs = {}
for i, fgroup in enumerate(fgroups):
    print(fgroup)
    with open('data/{}_R1_wbos.json'.format(fgroup), 'r') as f:
        wbos = json.load(f)

    wbos_filtered = []
    for wbo, name, smiles in wbos:
        if not 'phenoxide' in name:
            wbos_filtered.append([wbo, name, smiles])
    # Select molecules (up to 15) that are evenly distributed along WBOs
    wbos_dist = [round(wbo[0], 2) for wbo in wbos_filtered]
    differences = [j-i for i, j in zip(wbos_dist[:-1], wbos_dist[1:])]
    indices = [i for i, j in enumerate(differences) if j > 0]
    selected_mols = [wbos_filtered[i+1] for i in indices]
    selected_mols.append(wbos_filtered[-1])
    selected_mols.insert(0, wbos_filtered[0])


    if len(selected_mols) < 3:
        print('WBO range is too small. Dropping {}'.format(fgroup))
        continue
    if len(selected_mols) > 10:
        # remove some
        to_remove = selected_mols[1::int(round(len(selected_mols)/(len(selected_mols) - 10)))]
        wbo_filtered_2 = []
        for j in selected_mols:
            if j in to_remove:
                continue
            else:
                wbo_filtered_2.append(j)
        wbo_filtered_2.append(wbos_filtered[-1])

        #wbo_filtered = list(filter(lambda a: a in to_remove, wbos_filtered))
        selected_mols = wbo_filtered_2
    if selected_mols[-1][0] - selected_mols[-2][0] < 0.01:
        selected_mols.remove(selected_mols[-2])
    if selected_mols[1][0] - selected_mols[0][0] < 0.01:
        selected_mols.remove(selected_mols[1])

    smiles = [wbo[-1] for wbo in selected_mols]

    wbos = [wbo[0] for wbo in selected_mols]
    c_oe = openeye.oechem.OEColor(colors[color_keys[i]])
    fragmenter.chemi.to_pdf(smiles, fname='figures/qcarchive_torsiondrives/{}_to_drive.pdf'.format(fgroup), bo=wbos, bond_map_idx=(4, 5),
                 align=smiles[0], color=c_oe)


    print(fgroup, len(selected_mols))
    job_indices = []
    for wbo, name, smiles in selected_mols:
        mol =openeye.oechem.OEMol()
        openeye.oechem.OESmilesToMol(mol, smiles)
        mapped_mol = add_atom_map(mol)

        # Find bond
        central_bond = find_central_bond(mapped_mol)
        torsion = fragmenter.torsions.find_torsion_around_bond(mapped_mol, bond=central_bond)
        # get job index. Save job indices and WBOs for each fgroup - this will make analysis easier
        mapped_smiles = openeye.oechem.OEMolToSmiles(mapped_mol)
        job_index = cmiles.utils.to_canonical_label(openeye.oechem.OEMolToSmiles(mapped_mol), torsion)
        #job_index_smi = openeye.oechem.OEMolToSmiles(mapped_mol)
        job_index_smi=smiles
        if fgroup in to_replace:
            if isinstance(to_replace[fgroup], list):
                for index in to_replace[fgroup]:
                    if job_index == index[0]:
                        job_index = index[1]
            if job_index == to_replace[fgroup][0]:
                job_index = to_replace[fgroup][1]
        job_indices.append([job_index, torsion, name, wbo, job_index_smi])

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

