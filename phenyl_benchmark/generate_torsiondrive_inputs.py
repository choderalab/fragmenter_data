from openeye import oechem
import cmiles
import json
from fragmenter import torsions, chemi
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
    mol = oechem.OEMol(mol)
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


def generate_torsiondrive_input():
    """
    Generate input starting geometry, torsion to drive and some provenance
    Returns
    -------

    """
    pass


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

n = 0
for i, fgroup in enumerate(fgroups):
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
    if len(wbos_filtered) > 12:
        # remove some
        to_remove = wbos_filtered[2::int(len(wbos_filtered)/(len(wbos_filtered) - 12))]
        wbo_filtered_2 = []
        for j in wbos_filtered:
            if j in to_remove:
                continue
            else:
                wbo_filtered_2.append(j)

        #wbo_filtered = list(filter(lambda a: a in to_remove, wbos_filtered))
        wbos_filtered = wbo_filtered_2

    smiles = [wbo[-1] for wbo in wbos_filtered]

    wbos = [wbo[0] for wbo in wbos_filtered]
    c_oe = oechem.OEColor(colors[color_keys[i]])
    chemi.to_pdf(smiles, fname='figures/qcarchive_torsiondrives/{}_to_drive.pdf'.format(fgroup), bo=wbos, bond_map_idx=(4, 5),
                 align=smiles[0], color=c_oe)


    print(fgroup, len(wbos_filtered))
    n += len(wbos_filtered)
    td_input = {}
    job_indices = []
    for wbo, name, smiles in wbos_filtered:
        mol = oechem.OEMol()
        oechem.OESmilesToMol(mol, smiles)
        mapped_mol = add_atom_map(mol)

        # Find bond
        central_bond = find_central_bond(mapped_mol)
        torsion = torsions.find_torsion_around_bond(mapped_mol, bond=central_bond)
        # get job index. Save job indices and WBOs for each fgroup - this will make analysis easier
        mapped_smiles = oechem.OEMolToSmiles(mapped_mol)
        job_index = cmiles.utils.to_canonical_label(oechem.OEMolToSmiles(mapped_mol), torsion)
        job_indices.append([job_index, torsion, name, wbo])

    with open('data/{}_td_job_indices.json'.format(fgroup), 'w') as f:
        json.dump(job_indices, f, indent=2, sort_keys=True)

        # Prepare inputs in one large JSON

print(n)

