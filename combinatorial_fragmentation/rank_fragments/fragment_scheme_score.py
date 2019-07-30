import fragmenter
import json
from openeye import oechem
import matplotlib.pyplot as plt


def fragment_bond(mapped_mol, bond, threshold, path, functional_groups, keep_non_rotor):
    """
    Generate fragment for bond
    Parameters
    ----------
    mapped_mol :
    bond :

    Returns
    -------

    """
    print(functional_groups)
    print(bool(functional_groups))
    f = fragmenter.fragment.WBOFragmenter(mapped_mol, functional_groups)
    # Capture options
    f._options['keep_non_rotor_ring_substituents'] = keep_non_rotor

    # Add threshold as attribute because it is used in more than one function
    setattr(f, 'threshold', threshold)
    f._options.update({'threshold': threshold, 'path': path})
    # Calculate WBO for molecule
    #f.calculate_wbo()
    f._get_rotor_wbo()
    # Find ring systems
    f._find_ring_systems(keep_non_rotor_ring_substituents=keep_non_rotor)

    # Build fragment
    if bond not in f.rotors_wbo:
        bond = tuple(reversed(bond))
    f.build_fragment(bond, heuristic=path)

    return f

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="calculate OE WBO")
    parser.add_argument('-n', '--name', type=str, help='Molecule name with number of its state')
    parser.add_argument('-t', '--threshold', type=str, help='Threshold for fragmenter')
    parser.add_argument('-p', '--path', type=str, help="path for growing fragment")
    parser.add_argument('-f', '--functional_group', type=str)
    parser.add_argument('-nr', '--keep_non_rotor', type=str)
    args = parser.parse_args()
    name = args.name
    threshold = args.threshold
    path = args.path
    functional_groups = args.functional_group
    print(functional_groups)
    keep_non_rotor = args.keep_non_rotor
    if functional_groups == 'False':
        functional_groups = False
    if functional_groups == 'None':
        functional_groups = None
    if keep_non_rotor == 'False':
        keep_non_rotor = False
    if keep_non_rotor == 'True':
        keep_non_rotor = True
    threshold = float(threshold)

    with open('selected/{}/{}_selected_bonds.json'.format(name, name), 'r') as f:
        content = json.load(f)
    bonds = content['bonds']
    mapped_smiles = content['parent_smiles']

    with open('selected/{}/{}_frag_with_scores.json'.format(name, name), 'r') as f:
        frag_scores = json.load(f)
    with open('selected/{}/rescore/{}_frag_with_scores.json'.format(name, name), 'r') as f:
        frag_scores_2 = json.load(f)

    failures = {}
    mapped_mol = oechem.OEMol()
    oechem.OESmilesToMol(mapped_mol, mapped_smiles)
    charged_mol = fragmenter.chemi.get_charges(mapped_mol)
    score_size = {}
    for bond in bonds:
        ser_bond = fragmenter.utils.serialize_bond(bond)
        f = fragment_bond(charged_mol, tuple(bond), threshold, path, functional_groups, keep_non_rotor)
        frag_dict = f.to_json()
        frag_key = list(frag_dict.keys())[0]
        frag = frag_dict[frag_key]['cmiles_identifiers']['canonical_isomeric_smiles']
        frags = frag_scores[ser_bond]['frags']
        mmd_scores = frag_scores[ser_bond]['mmd_scores']
        if not frag in frags:
            frag = frag_key
            if not frag_key in frags:
                print('{} not found in bond {}'.format(frag, bond))
                failures[ser_bond] = frag
                continue

        idx = frags.index(frag)

        norm = plt.Normalize(min(mmd_scores), max(mmd_scores))
        normed_scores = norm(mmd_scores)
        score = mmd_scores[idx]
        normed_score = normed_scores[idx]

        frags_2 = frag_scores_2[ser_bond]['frags']
        mmd_scores_2 = frag_scores_2[ser_bond]['mmd_scores']
        idx_2 = frags_2.index(frag)

        norm_2 = plt.Normalize(min(mmd_scores_2), max(mmd_scores_2))
        normed_scores_2 = norm(mmd_scores_2)
        normed_score_2 = normed_scores_2[idx]

        if tuple(bond) not in f.fragments:
            bond = tuple(reversed(bond))
        mol = f.fragments[tuple(bond)]
        size = oechem.OECount(mol, oechem.OEIsHeavy())
        score_size[ser_bond] = [frag, score, normed_score, normed_score_2, size]
        score_size['provenance'] = frag_dict[frag_key]['provenance']
    filename = 'selected/{}/{}_{}_{}_{}_{}_score.json'.format(name, name, threshold,
                                                             path, functional_groups,
                                                             keep_non_rotor)

    with open(filename,'w') as f:
        json.dump(score_size, f, indent=2, sort_keys=True)

    if len(failures) > 0:

        filename = 'selected/{}/{}_{}_{}_{}_{}_failure.json'.format(name, name, threshold,
                                                                 path, functional_groups,
                                                                 keep_non_rotor)

        with open(filename,'w') as f:
            json.dump(failures, f, indent=2, sort_keys=True)

