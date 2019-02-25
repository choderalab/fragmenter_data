
# coding: utf-8

# In[175]:


#get_ipython().run_line_magic('matplotlib', 'inline')
import fragmenter
import json
import cmiles
import glob


# In[2]:
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="calculate OE WBO")
    parser.add_argument('-i', '--infile', type=str, help='Input JSON fragments file')
    parser.add_argument('-n', '--name', type=str, help='Molecule name')
    args = parser.parse_args()
    infile = args.infile
    name = args.name

#frag_jsons = glob.glob('../fragment/*fragments.json')
#frag_jsons.remove('../fragment/Trametinib_fragments.json')
    frag_jsons = [infile]

    # In[3]:
    with open(infile, 'r') as f:
        frags = json.load(f)

    # frags = {}
    # for f in frag_jsons:
    #     ki = f.split('/')[-1].split('_')[0]
    #     with open(f, 'r') as f:
    #         frags[ki] = json.load(f)


    # In[4]:


    omega_failures = []
    input_molecules = {}
    n = 0
    #for ki in frags:
        #print(ki)
    for frag in frags:
        # mapped_smiles
        qcarchive_input = {'type': 'optimization_input'}
        mol_id = frags[frag]['identifiers']
        provenance = frags[frag]['provenance']
        mapped_smiles = mol_id['canonical_isomeric_explicit_hydrogen_mapped_smiles']
        mapped_mol = cmiles.utils.load_molecule(mapped_smiles, toolkit='openeye')
        try:
            confs = fragmenter.chemi.generate_conformers(mapped_mol, dense=False)
        except:
            print('failed')
            print(frag)
            omega_failures.append(frag)
            # if name not in omega_failures:
            #     omega_failures[name] = []
            # omega_failures[name].append(frag)
            continue
        # Generate list of molecules for QC input
        qcschema_molecules = [cmiles.utils.mol_to_map_ordered_qcschema(conf, mol_id) for conf in confs.GetConfs()]
        print(len(qcschema_molecules))
        n += len(qcschema_molecules)
        qcarchive_input['initial_molecule'] = qcschema_molecules
        input_molecules[frag]= qcarchive_input
        input_molecules[frag]['provenance'] = frags[frag]['provenance']
    fname = 'mini_drug_bank/{}_bo_input.json'.format(name)
    with open(fname, 'w') as f:
        json.dump(input_molecules, f, indent=2, sort_keys=True)


# In[6]:


    print(omega_failures)

