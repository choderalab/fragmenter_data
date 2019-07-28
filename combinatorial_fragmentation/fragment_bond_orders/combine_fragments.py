"""
Combine fragments and cluster them by bonds to prep for scoring
"""

import json
import glob
from fragmenter import chemi, utils


def deserialize_dict(json_dict):
    deserialize_dict = {}
    for frag in json_dict:
        deserialize_dict[frag] = {}
        if 'map' in frag:
            deserialize_dict[frag] = json_dict[frag]
        else:
            for key in json_dict[frag]:
                if 'map' in key:
                    deserialize_dict[frag][key] = json_dict[frag][key]
                else:
                    des_key = utils.deserialize_bond(key)
                    deserialize_dict[frag][des_key] = json_dict[frag][key]
    return deserialize_dict

def check_for_key(bond_dict, key):
    if not key in bond_dict:
        # Try reverse
        rev_key = tuple(reversed(key))
        if not rev_key in bond_dict:
            return None
        else:
            return bond_dict[rev_key]
    else:
        return bond_dict[key]

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="calculate OE WBO")
    parser.add_argument('-n', '--name', type=str, help='Molecule name')
    args = parser.parse_args()
    name = args.name

    parent_file = 'validation_set/{}/{}_parent_oe_wbo.json'.format(name, name)
    all_frag_files = glob.glob('validation_set/{}/{}_frag*'.format(name, name))

    with open(parent_file, 'r') as f:
        parent_wbos = json.load(f)
    parent_wbos_des = deserialize_dict(parent_wbos)
    mapped_parent = parent_wbos_des['map_to_parent']
    mapped_parent_mol = chemi.smiles_to_oemol(mapped_parent, normalize=False)
    parent_key = [key for key in list(parent_wbos.keys()) if 'map' not in key][0]

    wbo_by_bond = {}
    for bond in mapped_parent_mol.GetBonds():
        if bond.IsInRing():
            continue
         # keep bonds that do not include H
        a1 = bond.GetBgn()
        a2 = bond.GetEnd()
        if a1.IsHydrogen() or a2.IsHydrogen():
            continue
        if a1.GetDegree() == 1 or a2.GetDegree() == 1:
            # Terminal
            continue
        bond_key = (a1.GetMapIdx(), a2.GetMapIdx())
        wbo_by_bond[bond_key] = {}
        try:
            wbo_by_bond[bond_key][parent_key] = {'individual_confs': parent_wbos_des[parent_key][bond_key]['individual_confs'],
                                                'elf_estimate': parent_wbos_des[parent_key][bond_key]['elf_estimate'],
                                                'map_to_parent': mapped_parent}
        except KeyError:
            rev_bond_key = tuple(reversed(bond_key))
            wbo_by_bond[bond_key][parent_key] = {'individual_confs': parent_wbos_des[parent_key][rev_bond_key]['individual_confs'],
                                                'elf_estimate': parent_wbos_des[parent_key][rev_bond_key]['elf_estimate'],
                                                'map_to_parent': mapped_parent}

    # Organize by bond
    for file in all_frag_files:
        with open(file, 'r') as f:
            frag = json.load(f)
        des_frag = deserialize_dict(frag)
        frag_key = list(des_frag.keys())[0]
        for bond in wbo_by_bond:
            wbos = check_for_key(des_frag[frag_key], bond)
            if not wbos:
                continue
            if frag_key not in wbo_by_bond[bond]:
                wbo_by_bond[bond][frag_key] = wbos
                wbo_by_bond[bond][frag_key]['map_to_parent'] = des_frag[frag_key]['"map_to_parent"']
            else:
                if des_frag[frag_key]['"map_to_parent"'] != wbo_by_bond[bond][frag_key]['map_to_parent']:
                    # Add this fragment. First check index
                    index = (file.split('/')[-1].split('_')[4])
                    if index[-1].isdigit():
                        frag_key_new = frag_key + '_' + index
                    else:
                        frag_key_new  = frag_key + '_1'
                    wbo_by_bond[bond][frag_key_new] = wbos
                    wbo_by_bond[bond][frag_key_new]['map_to_parent'] = des_frag[frag_key]['"map_to_parent"']
    # Serialize
    wbo_by_bond_ser = {}
    for bond in wbo_by_bond:
        ser_bond = utils.serialize_bond(bond)
        wbo_by_bond_ser[ser_bond] = wbo_by_bond[bond]
    with open('validation_set/{}/{}_oe_wbo_by_bond.json'.format(name, name), 'w') as f:
        json.dump(wbo_by_bond_ser, f, indent=2, sort_keys=True)



