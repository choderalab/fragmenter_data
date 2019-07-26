import fragmenter
from openeye import oechem
import pandas as pd
import json
import socket
import getpass

validation_set_input = pd.read_csv('../filter/validation_set/drugbank_filtered.csv')

# Clean up names = remove space because it causes bash scripts to crash
for name in validation_set_input['name']:
    if ' ' in name:
        name_split = name.split(' ')
        new_name = ''
        for n in name_split[:-1]:
            new_name += n 
            new_name += '_'
        new_name += name_split[-1]
        validation_set_input.replace(to_replace=name, value=new_name, inplace=True)
    if name[0].isdigit():
        new_name = 'four' + name[1:]
        name = new_name
        validation_set_input.replace(to_replace=name, value=new_name, inplace=True)
# save csv with fixed names
validation_set_input.to_csv('../filter/validation_set/drugbank_filtered_fixed_names.csv')

provenance = {}
provenance['creator'] = "fragmenter"
provenance['hostname'] = socket.gethostname()
provenance['username'] = getpass.getuser()
provenance['routine'] = {'enumerate_states' : {'version': fragmenter.__version__,
                        'options': {'stereoisomers': False, 'explicit_h': False, 'filter_nitro': True}}}

all_states =[]
all_names = []
for smiles, name in zip(validation_set_input.smiles, validation_set_input.name):
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, smiles)
    mol.SetTitle(name)
    
    states, names = fragmenter.fragment.enumerate_states(mol,stereoisomers=False, explicit_h=False, return_names=True,
                                                   filter_nitro=True)

    states_dict = {'provenance': provenance, 'states': states}
    states_dict['provenance']['routine']['enumerate_states']['parent_molecule_name'] = name
    states_dict['provenance']['routine']['enumerate_states']['parent_molecule'] = smiles
    with open('validation_set/{}_states.json'.format(name), 'w') as f:
        json.dump(states_dict, f, indent=2, sort_keys=True)
    all_states.extend(states)
    all_names.extend(names)

fragmenter.chemi.to_pdf(molecules=all_states, names=all_names, fname='enumerated_states.pdf')
print('Number of molelcules to fragment: {}'.format(len(all_states)))