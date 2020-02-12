"""
Fragment selected set using the Pfizer scheme from https://pubs.acs.org/doi/10.1021/acs.jcim.9b00373
For each fragment, calculate WBO distributions and save
"""

import fragmenter
import openeye
from openeye import oechem, oequacpac
import cmiles
import json
import itertools
import os

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
    with open('{}/{}_wbo_dists_fixed_1.json'.format(name, name), 'r') as f:
        already_calculated = json.load(f)

    # deserialize
    wbo_dists_des = {}
    for b in wbo_dists:
        b_des = fragmenter.utils.deserialize_bond(b)
        wbo_dists_des[b_des] = wbo_dists[b]

    parent_smiles = selected['parent_smiles']
    parent_mol = oechem.OEMol()
    oechem.OESmilesToMol(parent_mol, parent_smiles)

    fragmenter_engine = fragmenter.fragment.PfizerFragmenter(parent_mol)
    fragmenter_engine.fragment()
    fragmenter_engine.depict_fragments(fname='{}/{}_pfizer_frags.pdf'.format(name, name))

    # Get parent wbo dist
    frags = {}
    already_seen = {}
    for bond in selected['bonds']:
        already_seen[tuple(bond)] = {}
    for bond in selected['bonds']:
        t_bond = tuple(bond)
        if not t_bond in frags:
            frags[t_bond] = {}
        if not t_bond in wbo_dists_des:
            wbos = wbo_dists_des[(tuple(reversed(t_bond)))]
        else:
            wbos = wbo_dists_des[t_bond]
        for key in wbo_dists_des[t_bond]:
            if 'parent' in key:
                frags[t_bond]['parent'] = {'wbo_dist': wbo_dists_des[t_bond][key]['individual_confs'],
                                           'frag': wbo_dists_des[t_bond][key]['map_to_parent'],
                                           'elf10_wbo': wbo_dists_des[t_bond][key]['elf_estimate']}
        if not t_bond in fragmenter_engine.fragments:
            pfizer_frag = fragmenter_engine.fragments[(tuple(reversed(t_bond)))]
        else:
            pfizer_frag = fragmenter_engine.fragments[t_bond]

        frag = fragmenter.chemi.get_charges(pfizer_frag, strict_types=False, strict_stereo=False, keep_confs=-1)
        bo = get_bond(frag, t_bond)
        elf10_wbo = bo.GetData('WibergBondOrder')
        frags[t_bond] = {'frag': oechem.OEMolToSmiles(frag), 'wbo_dist': [], 'elf10_wbo': elf10_wbo}

        # Check if this fragment with this bond already has WBO calculated
        mol_copy = oechem.OEMol(pfizer_frag)
        cmiles.utils.remove_atom_map(mol_copy)
        smiles = oechem.OEMolToSmiles(mol_copy)
        if smiles in wbo_dists_des[t_bond]:
            frags[t_bond]['wbo_dist'] = wbo_dists_des[t_bond][smiles]['individual_confs']
        elif smiles in already_seen[t_bond]:
            print('already calculated {}'.format(smiles))
            frags[t_bond]['wbo_dist'] = already_seen[t_bond][smiles]
        else:
            print('{} not in dataset'.format(smiles))
            for b in already_seen:
                if smiles not in already_seen[b]:
                    already_seen[b][smiles] = []
            for i, conf in enumerate(frag.GetConfs()):
                print('{} out of {}'.format(i, frag.GetMaxConfIdx()))
                mol_copy = oechem.OEMol(conf)
                if oequacpac.OEAssignPartialCharges(mol_copy, oequacpac.OECharges_AM1BCCSym):
                    for b in already_seen:
                        bo = get_bond(mol_copy, b)
                        if not bo:
                            continue
                        else:
                            already_seen[b][smiles].append(bo.GetData('WibergBondOrder'))
            frags[t_bond]['wbo_dist'] = already_seen[t_bond][smiles]

    # serialize
    frags_ser = {}
    for bond in frags:
        b_ser = fragmenter.utils.serialize_bond(bond)
        frags_ser[b_ser] = frags[bond]

    frags_ser['provenance'] = {'fragmenter_version': fragmenter.__version__,
                               'openeye_version': openeye.__version__}
    try:
        os.mkdir('{}'.format(name))
    except FileExistsError:
        print('{} already exists. Files will be overwritten'.format(name))
    with open('{}/{}_pfizer_wbo_dists.json'.format(name, name), 'w') as f:
        json.dump(frags_ser, f, indent=2, sort_keys=True)


