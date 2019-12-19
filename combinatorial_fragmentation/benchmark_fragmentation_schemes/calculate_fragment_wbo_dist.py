"""
Calculate WBO distribution of fragments generated with different parameters of WBOFragmenter on selected
molecules and bonds. If the fragment WBO distribution was already calculated in the combinatorial fragmentation,
use that. Otherwise, generate conformers and calculate WBO distributions
"""

import fragmenter
from openeye import oechem, oequacpac
import cmiles
import json
import itertools
import os


def fragment_bond(mapped_mol, bond, threshold, path, functional_groups, keep_non_rotor, keep_confs=-1):
    """
    Fragment molecule around specified bond using provided parameters. Because it is fragmenting for one
    bond, this function mostly uses hidden function because it is not the canonical way to use fragmenter

    Parameters
    ----------
    mapped_mol : oemol with map indices
    bond : tuple of bond map indices
    threshold : float
    path : str
        options: wbo, path_length
    functional_groups : bool
        If False, do not keep functional groups. If None, use fragmenter defined functional groups
    keep_non_rotor : bool
        If True, keep non rotor substituents, If False, do not
    kwargs :  kwargs for get_charges

    Returns
    -------
    WBOFragmenter

    """
    f = fragmenter.fragment.WBOFragmenter(mapped_mol, functional_groups)
    # Capture options
    f._options['keep_non_rotor_ring_substituents'] = keep_non_rotor

    # Add threshold as attribute because it is used in more than one function
    setattr(f, 'threshold', threshold)
    f._options.update({'threshold': threshold, 'path': path})
    f._get_rotor_wbo()
    # Find ring systems
    f._find_ring_systems(keep_non_rotor_ring_substituents=keep_non_rotor)

    # Build fragment
    if bond not in f.rotors_wbo:
        bond = tuple(reversed(bond))
    f._build_fragment(bond, heuristic=path, strict_types=False, keep_confs=keep_confs, strict_stereo=False)

    return f


def get_bond(mol, bond_tuple):
    """
    Get bond in molecule
    Parameters
    ----------
    mol : oemole with map indices
    bond_tuple : tuple with map indices of bond

    Returns
    -------
    bond if found, False otherwise
    """

    a1 = mol.GetAtom(oechem.OEHasMapIdx(bond_tuple[0]))
    a2 = mol.GetAtom(oechem.OEHasMapIdx(bond_tuple[1]))
    if not a1 or not a2:
        return False
    bond = mol.GetBond(a1, a2)
    if not bond:
        return False
    return bond


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="calculate OE WBO")
    parser.add_argument('-n', '--name', type=str, help='Molecule name with number of its state')
    args = parser.parse_args()
    name = args.name

    # Get selected bonds from rank fragmnets and load in precalculated WBO distribution to save compute time
    with open('../rank_fragments/selected/{}/{}_selected_bonds.json'.format(name, name), 'r') as f:
        selected = json.load(f)
    with open('../rank_fragments/selected/{}/{}_oe_wbo_with_score.json'.format(name, name), 'r') as f:
        wbo_dists = json.load(f)

    # deserialize
    wbo_dists_des = {}
    for b in wbo_dists:
        b_des = fragmenter.utils.deserialize_bond(b)
        wbo_dists_des[b_des] = wbo_dists[b]

    # Get parent wbo dist
    frags = {}
    for bond in selected['bonds']:
        t_bond = tuple(bond)
        if not t_bond in frags:
            frags[t_bond] = {}
        if not t_bond  in wbo_dists_des:
            wbos = wbo_dists_des[(tuple(reversed(t_bond)))]
        else:
            wbos = wbo_dists_des[t_bond]
        for key in wbo_dists_des[t_bond]:
            if 'parent' in key:
                frags[t_bond]['parent'] = {'wbo_dist': wbo_dists_des[t_bond][key]['individual_confs'],
                                           'frag': wbo_dists_des[t_bond][key]['map_to_parent']}

    parent_smiles = selected['parent_smiles']
    parent_mol = oechem.OEMol()
    oechem.OESmilesToMol(parent_mol, parent_smiles)
    # Calculate WBOs for parent so it does not get calculated every time for fragmneter
    parent_mol = fragmenter.chemi.get_charges(parent_mol, strict_types=False, keep_confs=-1)

    parameters = {'threshold': [0.1, 0.01, 0.001],
                  'path': ['path_length', 'wbo'],
                  'functional_groups': [None, False],
                  'keep_non_rotor': [True, False]}

    already_seen = {}
    for bond in selected['bonds']:
        already_seen[tuple(bond)] = {}
    for bond in selected['bonds']:
        print(bond)
        t_bond = tuple(bond)
        for t, p, f, r in itertools.product(parameters['threshold'], parameters['path'],
                                            parameters['functional_groups'], parameters['keep_non_rotor']):

            key = '{}_{}_{}_{}'.format(str(t), p, str(f), str(r))
            parent_mol_copy = oechem.OEMol(parent_mol)
            frag = fragment_bond(parent_mol_copy, t_bond, threshold=t, path=p, functional_groups=f, keep_non_rotor=r,
                                 keep_confs=-1)
            if not t_bond in frag.fragments:
                f = frag.fragments[tuple(reversed(t_bond))]
            else:
                f = frag.fragments[t_bond]
            frags[t_bond][key] = {'frag': oechem.OEMolToSmiles(f), 'wbo_dist': []}
            mol_copy = oechem.OEMol(f)
            cmiles.utils.remove_atom_map(mol_copy)
            smiles = oechem.OEMolToSmiles(mol_copy)
            if smiles in wbo_dists_des[t_bond]:
                frags[t_bond][key]['wbo_dist'] = wbo_dists_des[t_bond][smiles]['individual_confs']
            elif smiles in already_seen[t_bond]:
                print('already calculated {}'.format(smiles))
                frags[t_bond][key]['wbo_dist'] = already_seen[t_bond][smiles]
            else:
                print('{} not in dataset'.format(smiles))
                for b in already_seen:
                    if smiles not in already_seen[b]:
                        already_seen[b][smiles] = []
                for conf in f.GetConfs():
                    mol_copy = oechem.OEMol(conf)
                    if oequacpac.OEAssignPartialCharges(mol_copy, oequacpac.OECharges_AM1BCCSym):
                        for b in already_seen:
                            bo = get_bond(mol_copy, b)
                            if not bo:
                                continue
                            else:
                                already_seen[b][smiles].append(bo.GetData('WibergBondOrder'))
                frags[t_bond][key]['wbo_dist'] = already_seen[t_bond][smiles]
    #serialize
    frags_ser = {}
    for bond in frags:
        b_ser = fragmenter.utils.serialize_bond(bond)
        frags_ser[b_ser] = frags[bond]
    os.mkdir('{}'.format(name))
    with open('{}/{}_wbos_dists.json'.format(name, name), 'w') as f:
        json.dump(frags_ser, f, indent=2, sort_keys=True)

