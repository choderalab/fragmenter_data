
# coding: utf-8

# In[175]:


get_ipython().run_line_magic('matplotlib', 'inline')
import fragmenter
import json
import cmiles
import glob


# In[2]:


frag_jsons = glob.glob('../fragment/*fragments.json')
frag_jsons.remove('../fragment/Trametinib_fragments.json')


# In[3]:


frags = {}
for f in frag_jsons:
    ki = f.split('/')[-1].split('_')[0]
    with open(f, 'r') as f:
        frags[ki] = json.load(f)


# In[4]:


omega_failures = {}
input_molecules = {}
n = 0
for ki in frags:
    print(ki)
    input_molecules[ki] = {}
    for frag in frags[ki]:  
        # mapped_smiles
        qcarchive_input = {'type': 'optimization_input'}
        mol_id = frags[ki][frag]['identifiers']
        provenance = frags[ki][frag]['provenance']
        mapped_smiles = mol_id['canonical_isomeric_explicit_hydrogen_mapped_smiles']
        mapped_mol = cmiles.utils.load_molecule(mapped_smiles, toolkit='openeye')
        try:
            confs = fragmenter.chemi.generate_conformers(mapped_mol, dense=False)
        except:
            print('failed')
            print(frag)
            if ki not in omega_failures:
                omega_failures[ki] = []
            omega_failures[ki].append(frag)
            continue
        # Generate list of molecules for QC input
        qcschema_molecules = [cmiles.utils.mol_to_map_ordered_qcschema(conf, mol_id) for conf in confs.GetConfs()]
        print(len(qcschema_molecules))
        n += len(qcschema_molecules)
        qcarchive_input['initial_molecule'] = qcschema_molecules
        input_molecules[ki][frag]= qcarchive_input
        input_molecules[ki][frag]['provenance'] = frags[ki][frag]['provenance']
    fname = '{}_bo_input.json'.format(ki)
    with open(fname, 'w') as f:
        json.dump(input_molecules[ki], f, indent=2, sort_keys=True)


# In[6]:


print(omega_failures)

