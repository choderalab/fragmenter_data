from fragmenter import chemi, fragment
import json
from cmiles.utils import mol_to_smiles, has_atom_map
from cmiles import get_molecule_ids
from openeye import oechem
import warnings
import copy

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

    Returns
    -------
    json_dict: dict
        dictionary containing provenance and fragments.

    """
    routine = 'enumerate_fragments'
    provenance = _get_provenance(workflow_id='validation_set', routine=routine)
    with open('../../validation_set_workflow.json', 'r') as f:
        options = json.load(f)['validation_set']['fragmenter']['enumerate_fragments']
    scheme = options['scheme']
    if 'functional_groups' in options:
        functional_groups = options['functional_groups']
    else:
        functional_groups = None
    options = options['options']

    parent_molecule = chemi.standardize_molecule(molecule, title)
    parent_molecule_smiles = mol_to_smiles(parent_molecule, isomeric=True, explicit_hydrogen=False,
                                                    mapped=False)
    provenance['routine']['enumerate_fragments']['parent_molecule_name'] = parent_molecule.GetTitle()
    provenance['routine']['enumerate_fragments']['parent_molecule'] = parent_molecule_smiles


    if scheme == 'combinatorial':
        fragment_engine = fragment.CombinatorialFragmenter(parent_molecule, functional_groups=functional_groups)
    elif scheme == 'wiberg_bond_order':
        fragment_engine = fragment.WBOFragmenter(parent_molecule, functional_groups=functional_groups)
    else:
        raise ValueError("Only combinatorial and wiberg_bond_order are supported fragmenting schemes")

    fragment_engine.fragment(**options)

    if generate_vis:
        fname = '{}.pdf'.format(parent_molecule.GetTitle())
        fragment_engine.depict_fragments(fname=fname)

    if not fragment_engine.fragments:
        warnings.warn("No fragments were generated for {}".format(parent_molecule_smiles))
        # ToDo: if combinatorial does not have fragments it means that there was no point to cut and if
        # wbo does not have fragments it means no internal rotatable bond was found but it might still have torsions
        # to drive. Not yet sure how to handle these cases
        return

    if mol_provenance:
        provenance['routine']['enumerate_states'] = mol_provenance['routine']['enumerate_states']

    # Generate identifiers for fragments
    fragments_json_dict = {}
    for frag in fragment_engine.fragments:
        mols = fragment_engine.fragments[frag]
        if not isinstance(mols, list):
            mols = [mols]
        # Currently cmiles only takes SMILES or JSON as input. But maybe we should allow oemols?
        smiles = mol_to_smiles(mols[0], mapped=False)
        identifiers = get_molecule_ids(smiles, canonicalization='openeye')
        frag_smiles = identifiers['canonical_isomeric_smiles']
        if frag_smiles not in fragments_json_dict:
            fragments_json_dict[frag_smiles] = {'identifiers': identifiers}
            fragments_json_dict[frag_smiles]['provenance'] = copy.deepcopy(provenance)
            fragments_json_dict[frag_smiles]['provenance']['canonicalization'] = identifiers.pop('provenance')
            for i, mol in enumerate(mols):
                if has_atom_map(mol):
                    parent_map_smiles = oechem.OEMolToSmiles(mol)
                    if i == 0:
                        fragments_json_dict[frag_smiles]['provenance']['routine']['enumerate_fragments']['map_to_parent'] \
                            = [parent_map_smiles]
                        if scheme == 'wiberg_bond_order':
                            fragments_json_dict[frag_smiles]['provenance']['routine']['enumerate_fragments']['central_rot_bond'] \
                                = [frag]
                    if i > 0:
                        fragments_json_dict[frag_smiles]['provenance']['routine']['enumerate_fragments']['map_to_parent'].append(parent_map_smiles)
        else:
            # This is an equivalent fragment from a different part of the parent molecule
            for mol in mols:
                if has_atom_map(mol):
                    # This is the map to the parent. Meanwhile I'm storing this information in provenance because it's
                    # helpful for analysis. This info is probably extraneous for regular use of fragmenter
                    parent_map_smiles = oechem.OEMolToSmiles(mol)
                    fragments_json_dict[frag_smiles]['provenance']['routine']['enumerate_fragments']['map_to_parent'].append(parent_map_smiles)
                    # Also store rotatable bond this is from
                    if scheme == 'wiberg_bond_order':
                        fragments_json_dict[frag_smiles]['provenance']['routine']['enumerate_fragments']['central_rot_bond'].append(frag)
    if json_filename:
        with open(json_filename, 'w') as f:
            json.dump(fragments_json_dict, f, indent=2, sort_keys=True)

    return fragments_json_dict

def _get_provenance(workflow_id, routine):
    """
    Get provenance with keywords for routine

    Parameters
    ----------
    routine: str
        routine to get provenance for. Options are 'enumerate_states', 'enumerate_fragments', and 'generate_crank_jobs'
    options: str, optional. Default is None
        path to yaml file containing user specified options.

    Returns
    -------
    provenance: dict
        dictionary with provenance and routine keywords.

    """
    import fragmenter
    import uuid
    import socket
    import getpass
    fragmenter_version = fragmenter.__version__
    provenance = {'creator': fragmenter.__package__,
                  'job_id': str(uuid.uuid4()),
                  'hostname': socket.gethostname(),
                  'username': getpass.getuser(),
                  'workflow_id': workflow_id,
                  'routine': {routine: {
                      'version': fragmenter_version
                  }}}
    # # check version against version in wf
    # version_in_wf = self.off_workflow.get_options(routine)['version']
    # if version_in_wf != fragmenter_version:
    #     warnings.warn('You are using version {} of fragmenter for {}. The version specified in the registered QCArchive '
    #                   'workflow is {} '.format(fragmenter_version, routine, version_in_wf))

    return provenance


# In[2]:
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="fragment validation set")
    parser.add_argument('-i', '--infile', type=str, help='Input JSON fragments file')
    parser.add_argument('-n', '--name', type=str, help='Molecule name')
    args = parser.parse_args()
    infile = args.infile
    name = args.name

#frag_jsons = glob.glob('../fragment/*fragments.json')
#frag_jsons.remove('../fragment/Trametinib_fragments.json')
    #frag_jsons = [infile]

    # In[3]:
    with open(infile, 'r') as f:
        states = json.load(f)

    for i, state in enumerate(states['states']):
        title = '{}_{}'.format(name, str(i))
        json_filename = 'fragments_{}.json'.format(title)
        enumerate_fragments(state, title=title, mol_provenance=states['provenance'],
                                               json_filename=json_filename, generate_vis=True)