from fragmenter import chemi, fragment
import json

def enumerate_fragments(molecule, title='', mol_provenance=None, json_filename=None,
                        generate_vis=True):
    """
    Fragment molecule

    Parameters
    ----------
    molecule: Input molecule. Very permissive. Can be anything that OpenEye can parse
        SMILES string of molecule to fragment
    title: str, optional. Default empty str
        The title or name of the molecule. If empty stirng will use the IUPAC name for molecule title.
    mol_provenance: dict, optional. Default is None
        provenance for molecule. If the molecule is a state from enumerate_states, the provenance from enumerate_states
        should be used
    json_filename: str, optional. Default None
        If a filename is provided, will write output to json file.
    generate_vis: bool, optional, default False
        If True, will generate visualization of fragments from parent molecule

    """
    parent_molecule = chemi.smiles_to_oemol(molecule, name=title)
    fragment_engine = fragment.CombinatorialFragmenter(parent_molecule, functional_groups=False)
    fragment_engine.fragment()

    frag_json = fragment_engine.to_json()
    with open(json_filename, 'w') as f:
        json.dump(frag_json, f, indent=2, sort_keys=True)

    if len(fragment_engine.new_stereo) > 0:
        # Save fragments that are missing stereo because of new stereo center
        chemi.smiles_to_smi(fragment_engine.new_stereo, filename='{}_new_stereo_missing.smi'.format(title))


    if generate_vis:
        fname = '{}.pdf'.format(parent_molecule.GetTitle())
        fragment_engine.depict_fragments(fname=fname)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="fragment validation set")
    parser.add_argument('-i', '--infile', type=str, help='Input JSON fragments file')
    parser.add_argument('-n', '--name', type=str, help='Molecule name')
    args = parser.parse_args()
    infile = args.infile
    name = args.name

    with open(infile, 'r') as f:
        states = json.load(f)

    for i, state in enumerate(states['states']):
        title = '{}_{}'.format(name, str(i))
        json_filename = '{}_fragments.json'.format(title)
        enumerate_fragments(state, title=title, mol_provenance=states['provenance'],
                                               json_filename=json_filename, generate_vis=True)