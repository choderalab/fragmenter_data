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
        if 'WibergBondOrder' not in bond.GetData():
            warnings.warn("Molecule is missing WBO. It wasn't charged.")
            return False
        wbos[bond_key] = bond.GetData('WibergBondOrder')
    return wbos

def compute_wbos(map_to_parent):
    oemol = cmiles.utils.load_molecule(map_to_parent)
    try:
        charged = fragmenter.chemi.get_charges(oemol, keep_confs=-1, strict_types=False)
    except:
        print('Cannot charge {}'.format(oechem.OEMolToSmiles(oemol)))
        return False
    # Collect all wbos
    elf_wbo_estimate = collect_wbos(charged)
    if not elf_wbo_estimate:
        warnings.warn('{} was not charged'.format(oechem.OEMolToSmiles(charged)))
        return False
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

def serialize_key(key):
        if isinstance(key, (int, float)):
            key = (int(key), )

        return json.dumps(key)
def serialize(oe_wbo):
    serialized_wbo = {}
    for frag in oe_wbo:
        serialized_wbo[frag] = {}
        for bond in oe_wbo[frag]:
            key = serialize_key(bond)
            serialized_wbo[frag][key] = oe_wbo[frag][bond]
    return serialized_wbo

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="calculate OE WBO")
    parser.add_argument('-i', '--infile', type=str, help='Input JSON conformers file')
    parser.add_argument('-n', '--name', type=str, help='Molecule name')
    args = parser.parse_args()
    infile = args.infile
    name = args.name

    with open(infile, 'r') as f:
        frags = json.load(f)

    frag = list(frags.keys())[0]
    mapped_parent_smiles = frags[frag]['provenance']['routine']['fragment_molecule']['mapped_parent_smiles']
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
    try:
        os.mkdir('../../fragment_conformers/validation_set/{}'.format(name))
    except FileExistsError:
        print('{} directory already exists. Files will be overwritten'.format(name))
    fname = '../../fragment_conformers/validation_set/{}/{}_parent_conformers.json'.format(name, name)

    with open(fname, 'w') as f:
        json.dump(qcarchive_input, f, indent=2, sort_keys=True)

    # Save wbos

    wbos = {'{}_parent'.format(cmiles_identifiers['canonical_isomeric_smiles']): parent_wbos}
    try:
        os.mkdir(name)
    except FileExistsError:
         print('{} directory already exists. Files will be overwritten'.format(name))
    wbos = serialize(wbos)
    wbos['map_to_parent'] = mapped_parent_smiles


    with open('{}/{}_parent_oe_wbo.json'.format(name, name), 'w') as f:
        json.dump(wbos, f, indent=2, sort_keys=True)

    # Get wbo and conformers for fragments
    with open(infile, 'r') as f:
        fragments = json.load(f)
    #all_wbos = {}
    #all_conformers = {}
    oe_failures = []

    for i, frag in enumerate(fragments):
        wbos = {}
        conformers = {}
        mol_id = fragments[frag]['cmiles_identifiers']
        provenance = fragments[frag]['provenance']
        map_to_parents = fragments[frag]['provenance']['routine']['fragment_molecule']['map_to_parent']
        for j,  map_to_parent in enumerate(map_to_parents):
            if i > 0:
                frag_name = '{}_{}'.format(frag, str(j))
            else:
                frag_name = frag
            map_to_parent = map_to_parents[0]
            computed = compute_wbos(map_to_parent)
            if not computed:
                warnings.warn('Failed to charge molecule {}'.format(frag))
                oe_failures.append(frag)
                continue
            charged = computed[0]
            wbos[frag] = computed[1]
            wbos[frag]['map_to_parent'] = map_to_parent
            mapped_smiles = mol_id['canonical_isomeric_explicit_hydrogen_mapped_smiles']
            qcschema_molecules = [cmiles.utils.mol_to_map_ordered_qcschema(conf, mapped_smiles) for conf in charged.GetConfs()]
            conformers[frag] = {'initial_molecules': qcschema_molecules,
                                         'cmiles_identifiers': mol_id,
                                         'provenance': provenance}
            fname = '{}/{}_frag_{}_{}_oe_wbo.json'.format(name, name, str(i), str(j))
            wbos = serialize(wbos)
            with open(fname, 'w') as f:
                json.dump(wbos, f, sort_keys=True, indent=2)
            # Save conformers
            fname = '../../fragment_conformers/validation_set/{}/{}_{}_conformers.json'.format(name, name, str(i))
            with open(fname, 'w') as f:
                json.dump(conformers, f, sort_keys=True, indent=2)


    # Save failures
    if len(oe_failures) > 0:
        fname = '{}/{}_oe_failed_fragments.json'.format(name, name)
        with open(fname, 'w') as f:
            json.dump(oe_failures, f, indent=2, sort_keys=True)
