import fragmenter
import json
import cmiles
from openeye import oechem, oequacpac
import os
import warnings


def collect_wbos(molecule):
    wbos = {}
    for bond in molecule.GetBonds():
        bond_key = (bond.GetBgn().GetMapIdx(), bond.GetEnd().GetMapIdx())
        wbos[bond_key] = bond.GetData('WibergBondOrder')
    return wbos

def compute_wbos(map_to_parent):
    oemol = cmiles.utils.load_molecule(map_to_parent)
    try:
        charged = fragmenter.chemi.get_charges(oemol, keep_confs=-1)
    except:
        print('Cannot charge {}'.format(oechem.OEMolToSmiles(charged)))
        return False
    # Collect all wbos
    elf_wbo_estimate = collect_wbos(charged)
    wbos = {key: {'elf_estimate': elf_wbo_estimate[key], 'individual_confs': []} for key in elf_wbo_estimate}
    # Compute WBO for each conformer
    for i, conf in enumerate(charged.GetConfs()):
        mol_copy = oechem.OEMol(conf)
        # Get WBO
        if oequacpac.OEAssignPartialCharges(mol_copy, oequacpac.OECharges_AM1BCCSym):
            ind_wbos = collect_wbos(molecule=mol_copy)
            for b in ind_wbos:
                if b not in wbos and 0 in b:
                    wbos[b] = {'individual_confs': [ind_wbos[b]]}
                elif b not in wbos:
                    reverse = tuple(reversed(b))
                    wbos[reverse]['individual_confs'].append(ind_wbos[b])

                wbos[b]['individual_confs'].append(ind_wbos[b])
        else:
            print('AM1BCC charging failed for {}, {}'.format(str(i), i))
    return charged, wbos

def serialize(oe_wbo):
    serialized_wbo = {}
    for frag in oe_wbo:
        serialized_wbo[frag] = {}
        for bond in oe_wbo[frag]:
            key = fragmenter.workflow_api.serialize_key(bond)
            serialized_wbo[frag][key] = oe_wbo[frag][bond]
    return serialized_wbo

def organize_wbo(mapped_parent_smiles, oe_wbo, fragment_inputs):

    mapped_parent_mol = cmiles.utils.load_molecule(mapped_parent_smiles)
    # Find all relevant bonds
    rot_bonds = []
    for bond in mapped_parent_mol.GetBonds():
        if bond.IsInRing():
            continue
        # keep bonds that do not include H
        a1 = bond.GetBgn()
        a2 = bond.GetEnd()
        if a1.IsHydrogen() or a2.IsHydrogen():
            continue
        key = (bond.GetBgn().GetMapIdx(), bond.GetEnd().GetMapIdx())
        rot_bonds.append(key)
    frag_with_bond = {b:{} for b in rot_bonds }

    for frag in oe_wbo:
        frag_name = frag.split('_')
        if len(frag_name) > 1:
            index = frag_name[-1]
        else:
            index = 0
        frag_name = frag_name[0]
        if not index == 'parent':
            map_to_parent = fragment_inputs[frag_name]['provenance']['routine']['enumerate_fragments']['map_to_parent'][int(index)]
            oemol = cmiles.utils.load_molecule(map_to_parent)
        elif index == 'parent':
            oemol = cmiles.utils.load_molecule(mapped_parent_smiles)
        #r_bonds = [(b.GetBgn().GetMapIdx(), b.GetEnd().GetMapIdx()) for b in oemol.GetBonds() if b.IsRotor()]
        r_bonds = []
        for bond in oemol.GetBonds():
            if bond.IsInRing():
                continue
            # keep bonds that do not include H
            a1 = bond.GetBgn()
            a2 = bond.GetEnd()
            if a1.IsHydrogen() or a2.IsHydrogen():
                continue
            key = (bond.GetBgn().GetMapIdx(), bond.GetEnd().GetMapIdx())
            r_bonds.append(key)
        for b in r_bonds:
            reverse = tuple(reversed(b))
            if b not in oe_wbo[frag]:
                # Try reverse
                bo = oe_wbo[frag][reverse]
            else:
                bo = oe_wbo[frag][b]
            try:
                frag_with_bond[b][frag] = bo
                frag_with_bond[b][frag]['map_to_parent'] = map_to_parent
            except KeyError:
                frag_with_bond[reverse][frag] = bo
                frag_with_bond[reverse][frag]['map_to_parent'] = map_to_parent
    serialized = {}
    for bond in frag_with_bond:
        key = fragmenter.workflow_api.serialize_key(bond)
        serialized[key] = frag_with_bond[bond]
    return serialized



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="calculate OE WBO")
    parser.add_argument('-i', '--infile', type=str, help='Input JSON conformers file')
    parser.add_argument('-n', '--name', type=str, help='Molecule name')
    args = parser.parse_args()
    infile = args.infile
    name = args.name

    # Get parent
    index = name.split('_')[-1]
    with open('../../enumerate_states/states_{}.json'.format(name.split('_')[0]), 'r') as f:
        states = json.load(f)['states']
    mapped_parent_smiles = states[int(index)]
    mapped_parent_mol = oechem.OEMol()
    oechem.OESmilesToMol(mapped_parent_mol, mapped_parent_smiles)
    cmiles_identifiers = cmiles.get_molecule_ids(mapped_parent_smiles)
    # Check that both mapped SMILES are the same
    if not mapped_parent_smiles == cmiles_identifiers['canonical_isomeric_explicit_hydrogen_mapped_smiles']:
        warnings.warn('mapped smiles do not match. Replacing with original mapped SMILES')
        cmiles_identifiers['canonical_isomeric_explicit_hydrogen_mapped_smiles'] = mapped_parent_smiles

    # compute parent WBO (and generate conformers)
    charged, parent_wbos = compute_wbos(mapped_parent_smiles)
    qcschema_molecules = [cmiles.utils.mol_to_map_ordered_qcschema(conf, mapped_parent_smiles) for conf in charged.GetConfs()]
    # Save parent conformers
    qcarchive_input = {'initial_molecules': qcschema_molecules,
                       'cmiles_identifiers': cmiles_identifiers}
    fname = '../../fragment_conformers/validation_set/{}_parent_conformers.json'.format(name)

    with open(fname, 'w') as f:
        json.dump(qcarchive_input, f, indent=2, sort_keys=True)

    # Get wbo and conformers for fragments
    with open(infile, 'r') as f:
        fragments = json.load(f)
    all_wbos = {}
    all_conformers = {}
    for frag in fragments:
        mol_id = fragments[frag]['identifiers']
        provenance = fragments[frag]['provenance']
        map_to_parents = fragments[frag]['provenance']['routine']['enumerate_fragments']['map_to_parent']
        for i, map_to_parent in enumerate(map_to_parents):
            if i > 0:
                frag_name = '{}_{}'.format(frag, str(i))
            else:
                frag_name = frag

            charged, all_wbos[frag] = compute_wbos(map_to_parent)
            all_wbos[frag]['map_to_parent'] = map_to_parent
            mapped_smiles = mol_id['canonical_isomeric_explicit_hydrogen_mapped_smiles']
            qcschema_molecules = [cmiles.utils.mol_to_map_ordered_qcschema(conf, mapped_smiles) for conf in charged.GetConfs()]
            all_conformers[frag_name] = {'initial_molecules': qcschema_molecules,
                                         'cmiles_identifiers': mol_id,
                                         'provenance': provenance}
    # Save conformers
    fname = '../../fragment_conformers/validation_set/{}_conformers.json'.format(name)
    with open(fname, 'w') as f:
        json.dump(all_conformers, f, sort_keys=True, indent=2)

    # Save all oe wbo
    all_wbos['{}_parent'.format(cmiles_identifiers['canonical_isomeric_smiles'])] = parent_wbos
    serialized_wbo = serialize(all_wbos)
    try:
        os.mkdir(name)
    except FileExistsError:
         print('{} directory already exists. Files will be overwritten'.format(name))
    with open('{}/{}_oe_wbo.json'.format(name, name), 'w') as f:
        json.dump(serialized_wbo, f, indent=2, sort_keys=True)

    # Organize wbos for easier analysis
    organized_wbo = organize_wbo(mapped_parent_smiles, all_wbos, fragments)
    fname = '{}/{}_oe_wbo_by_bond.json'.format(name, name)
    with open(fname, 'w') as f:
        json.dump(organized_wbo, f, indent=2, sort_keys=True)
