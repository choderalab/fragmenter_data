import fragmenter
import cmiles
import json
import os
from openeye import oechem, oequacpac
import openeye

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

    with open('../enumerate_states/validation_set/{}_states.json'.format(name), 'r') as f:
        states = json.load(f)['states']

    for i, state in enumerate(states):
        new_name = '{}_{}'.format(name, i)
        with open('../fragment_bond_orders/validation_set/{}/{}_oe_wbo_by_bond.json'.format(new_name, new_name), 'r') as f:
            combinatorial_results = json.load(f)
        # Deserialize
        combinatorial_results_des = {}
        for bond in combinatorial_results:
            bond_des = fragmenter.utils.deserialize_bond(bond)
            combinatorial_results_des[bond_des] = combinatorial_results[bond]
        wbo_dists = {}
        print(state)
        mol = fragmenter.chemi.smiles_to_oemol(state)
        frag = fragmenter.fragment.PfizerFragmenter(mol)
        frag.fragment()

        for bond in frag.fragments:
            # First check if wbos were already calculated in combinatorial fragmentation
            if not bond in combinatorial_results_des:
                results = combinatorial_results_des[tuple(reversed(bond))]
            else:
                results = combinatorial_results_des[bond]

            # Find parent
            for m in results:
                if 'parent' in m:
                    parent_elf10_wbo = results[m]['elf_estimate']
                    parent_wbo_dist = results[m]['individual_confs']
            print(results.keys())
            fragment = frag.fragments[bond]
            fragment_copy = oechem.OEMol(fragment)
            cmiles.utils.remove_atom_map(fragment_copy)
            smiles = oechem.OEMolToSmiles(fragment_copy)
            if bond not in wbo_dists:
                wbo_dists[bond] = {'wbo_dist': [], 'frag': oechem.OEMolToSmiles(fragment)}
            if smiles in results:
                elf10_wbo = results[smiles]['elf_estimate']
                wbo_dists[bond]['elf10_wbo'] = elf10_wbo
                wbo_dists[bond]['wbo_dist'] = results[smiles]['individual_confs']
            else:
                print("{} not found in results".format(smiles))
                charged = fragmenter.chemi.get_charges(frag.fragments[bond], strict_types=False, keep_confs=-1)
                print(bond)
                oe_bond = get_bond(charged, bond)
                print(oe_bond)
                elf10_wbo = oe_bond.GetData('WibergBondOrder')
                wbo_dists[bond]['elf10_wbo'] = elf10_wbo
                for conf in charged.GetConfs():
                    mol_copy = oechem.OEMol(conf)
                    if oequacpac.OEAssignPartialCharges(mol_copy, oequacpac.OECharges_AM1BCCSym):
                        bo = get_bond(mol_copy, bond)
                        wbo_dists[bond]['wbo_dist'].append(bo.GetData('WibergBondOrder'))
            wbo_dists[bond]['parent_elf10_wbo'] = parent_elf10_wbo
            wbo_dists[bond]['parent_wbo_dist'] = parent_wbo_dist

        # serialize
        wbo_dists_ser = {}
        for bond in wbo_dists:
            bond_ser = fragmenter.utils.serialize_bond(bond)
            wbo_dists_ser[bond_ser] = wbo_dists[bond]
        wbo_dists_ser['provenance'] = {'fragmenter_version': fragmenter.__version__, 'openeye_versoin': openeye.__version__}

        try:
            os.mkdir('{}'.format(new_name))
        except FileExistsError:
            print('{} already exists. Files will be overwritten'.format(new_name))
        with open('{}/{}_pfizer_fragment_wbo.json'.format(new_name, new_name), 'w') as f:
            json.dump(wbo_dists_ser, f, indent=2, sort_keys=True)
        frag.depict_fragments('{}/{}_fragments.pdf'.format(new_name, new_name))
