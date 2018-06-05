
# coding: utf-8

# In[1]:

get_ipython().magic('matplotlib inline')
import matplotlib.pyplot as plt
from openeye import oechem, oedepict
from fragmenter import utils
import json
import numpy as np
from openmoltools import openeye
import itertools
import glob
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.image as mpimg
import os


# In[2]:

# load all directories
directories = [x[0] for x in os.walk(os.getcwd())]


# In[19]:

all_data = {}
for kinase_inhibitor in directories[2:]:
    ki_key = kinase_inhibitor.split('/')[-1]
    all_data[ki_key] = []
    output_files = glob.glob('{}/*.output.json'.format(ki_key))
    for file in output_files:
        f = open(file, 'r')
        data = json.load(f)
        f.close
        # Check for errors
        error = data['error']
        if error:
            print(error)
            continue
        data['bond_orders'] = utils.bond_order_from_psi4_raw_output(data['raw_output'])
        data.pop('raw_output')
        all_data[ki_key].append(data)


# In[21]:

for ki in all_data: 
    print(ki)
    conformations = len(all_data[ki])
    molecule, atom_map = utils.get_atom_map(tagged_smiles=all_data[ki][0]['tagged_smiles'])
    molecule_charged = openeye.get_charges(molecule)
    n_atoms = molecule.NumAtoms()
    # create bond dictionary (Skip all hydrogen bonds)
    bond_order_wiberg = {}
    bond_order_mayer = {}
    for bond in molecule.GetBonds():
        atom_1 = bond.GetBgn()
        atom_2 = bond.GetEnd()
        if atom_1.IsHydrogen() or atom_2.IsHydrogen():
            continue
        map_1 = atom_1.GetMapIdx()
        map_2 = atom_2.GetMapIdx()
        bond_order_wiberg[(map_1, map_2)] = np.zeros(conformations)
        bond_order_mayer[(map_1, map_2)] = np.zeros(conformations)

    # Populate array
    for k, data in enumerate(all_data[ki]):
        wiberg = data['bond_orders']['Wiberg_psi4']
        mayer = data['bond_orders']['Mayer_psi4']
        for i, j in bond_order_wiberg:
            bond_order_wiberg[(i, j)][k] = wiberg[i-1][j-1]
            bond_order_mayer[(i, j)][k] = mayer[i-1][j-1]
        
    # Generate figure with labeled map
    utils.png_atoms_labeled(all_data[ki][0]['tagged_smiles'], '{}_mapped.png'.format(ki), width=400, height=400, 
                            label_scale=2.0, scale_bondwidth=True)
    #Charge molecule
    # plot bond orders
    with PdfPages('{}_bond_orders.pdf'.format(ki)) as pdf:
        n = len(bond_order_wiberg)
        plots = 8
        chunk = n/(plots-2)
        keys = list(bond_order_wiberg.keys())
        #for i in range(chunk):
        #    for j in range(6)
        extra = n % (plots-2)
        for i in range(int(chunk)):
            for j in range(plots-2):
                try:
                    key = keys[(plots-2)*i+j]
                except IndexError:
                    continue
                figure = plt.subplot2grid((4,2), (0, 0), rowspan=2)
                img = mpimg.imread('{}_mapped.png'.format(ki))
                imgplot = plt.imshow(img, interpolation='none')
                plt.xticks([])
                plt.yticks([])
                bond = bond_order_wiberg[key]
                atom_1 = molecule_charged.GetAtom(oechem.OEHasMapIdx(key[0]))
                atom_2 = molecule_charged.GetAtom(oechem.OEHasMapIdx(key[1]))
                bond_oe = molecule_charged.GetBond(atom_1, atom_2)
                wiberg_oe = bond_oe.GetData('WibergBondOrder')
                if j == 0:
                    plt.subplot2grid((4, 2), (0, 1))
                    plt.hist(bond_order_wiberg[key], alpha=0.5, label='Wiberg_psi4')
                    plt.hist(bond_order_mayer[key], alpha=0.5, label='Mayer_psi4')
                    plt.vlines(wiberg_oe, ymin=0, ymax=20, label='Wiberg_OE')
                    plt.title(key)
                    lgd = plt.legend(prop={'size': 8}, loc='center left', bbox_to_anchor=(1, 0.5))
                if j == 2:
                    plt.subplot(plots/2, 2, 3)
                if j == 1:
                    plt.subplot(plots/2, 2, 4)
                else:
                    plt.subplot(plots/2, 2, j+3)
                plt.hist(bond_order_wiberg[key], alpha=0.5, label='Wiberg_psi4')
                plt.hist(bond_order_mayer[key], alpha=0.5, label='Mayer_psi4')
                plt.vlines(wiberg_oe, ymin=0, ymax=20, label='Wiberg_OE')
                #plt.xlim(0.9, 2.0)
                plt.title(key)
                plt.tight_layout()
            pdf.savefig(dpi=300, bbox_extra_artist=(lgd,), bbox_inches='tight')
            plt.close()
        if extra:
            figure = plt.subplot2grid((4,2), (0, 0), rowspan=2)
            imgplot = plt.imshow(img, interpolation='none')
            plt.xticks([])
            plt.yticks([])
            bond = bond_order_wiberg[key]
            atom_1 = molecule_charged.GetAtom(oechem.OEHasMapIdx(key[0]))
            atom_2 = molecule_charged.GetAtom(oechem.OEHasMapIdx(key[1]))
            bond_oe = molecule_charged.GetBond(atom_1, atom_2)
            wiberg_oe = bond_oe.GetData('WibergBondOrder')
            for m in range(extra):
                key = keys[((plots-2)*i+j) + m+1]
                plt.subplot(plots/2, 2, 4+m)
                plt.hist(bond_order_wiberg[key], alpha=0.5, label='Wiberg_psi4')
                plt.hist(bond_order_mayer[key], alpha=0.5, label='Mayer_psi4')
                plt.vlines(wiberg_oe, ymin=0, ymax=20, label='Wiberg_OE')
                if m == 0:
                    lgd = plt.legend(prop={'size': 8}, loc='center left', bbox_to_anchor=(1, 0.5))
                #plt.xlim(0.9, 2.0)
                plt.title(key)
                plt.tight_layout()
            pdf.savefig(dpi=300, bbox_extra_artist=(lgd,), bbox_inches='tight')
            plt.close()



# In[39]:

for ki in all_data: 
    print(ki)
    conformations = len(all_data[ki])
    molecule, atom_map = utils.get_atom_map(tagged_smiles=all_data[ki][0]['tagged_smiles'])
    molecule_charged = openeye.get_charges(molecule)
    n_atoms = molecule.NumAtoms()
    # create bond dictionary (Skip all hydrogen bonds)
    bond_order_wiberg = {}
    bond_order_mayer = {}
    for bond in molecule.GetBonds():
        atom_1 = bond.GetBgn()
        atom_2 = bond.GetEnd()
        if atom_1.IsHydrogen() or atom_2.IsHydrogen():
            continue
        if atom_1.IsInRing() and atom_2.IsInRing():
            continue
        map_1 = atom_1.GetMapIdx()
        map_2 = atom_2.GetMapIdx()
        bond_order_wiberg[(map_1, map_2)] = np.zeros(conformations)
        bond_order_mayer[(map_1, map_2)] = np.zeros(conformations)

    # Populate array
    for k, data in enumerate(all_data[ki]):
        wiberg = data['bond_orders']['Wiberg_psi4']
        mayer = data['bond_orders']['Mayer_psi4']
        for i, j in bond_order_wiberg:
            bond_order_wiberg[(i, j)][k] = wiberg[i-1][j-1]
            bond_order_mayer[(i, j)][k] = mayer[i-1][j-1]
        
    # Generate figure with labeled map
    utils.png_atoms_labeled(all_data[ki][0]['tagged_smiles'], '{}_mapped.png'.format(ki), width=400, height=400, 
                            label_scale=2.0, scale_bondwidth=True)
    #Charge molecule
    # plot bond orders
    with PdfPages('{}_bond_orders_no_rings.pdf'.format(ki)) as pdf:
        n = len(bond_order_wiberg)
        plots = 8
        chunk = n/(plots-2)
        keys = list(bond_order_wiberg.keys())
        #for i in range(chunk):
        #    for j in range(6)
        extra = n % (plots-2)
        for i in range(int(chunk)):
            for j in range(plots-2):
                try:
                    key = keys[(plots-2)*i+j]
                except IndexError:
                    continue
                figure = plt.subplot2grid((4,2), (0, 0), rowspan=2)
                img = mpimg.imread('{}_mapped.png'.format(ki))
                imgplot = plt.imshow(img, interpolation='none')
                plt.xticks([])
                plt.yticks([])
                bond = bond_order_wiberg[key]
                atom_1 = molecule_charged.GetAtom(oechem.OEHasMapIdx(key[0]))
                atom_2 = molecule_charged.GetAtom(oechem.OEHasMapIdx(key[1]))
                bond_oe = molecule_charged.GetBond(atom_1, atom_2)
                wiberg_oe = bond_oe.GetData('WibergBondOrder')
                if j == 0:
                    plt.subplot2grid((4, 2), (0, 1))
                    plt.hist(bond_order_wiberg[key], alpha=0.5, label='Wiberg_psi4')
                    plt.hist(bond_order_mayer[key], alpha=0.5, label='Mayer_psi4')
                    plt.vlines(wiberg_oe, ymin=0, ymax=20, label='Wiberg_OE')
                    plt.title(key)
                    lgd = plt.legend(prop={'size': 8}, loc='center left', bbox_to_anchor=(1, 0.5))
                if j == 2:
                    plt.subplot(plots/2, 2, 3)
                if j == 1:
                    plt.subplot(plots/2, 2, 4)
                else:
                    plt.subplot(plots/2, 2, j+3)
                plt.hist(bond_order_wiberg[key], alpha=0.5, label='Wiberg_psi4')
                plt.hist(bond_order_mayer[key], alpha=0.5, label='Mayer_psi4')
                plt.vlines(wiberg_oe, ymin=0, ymax=20, label='Wiberg_OE')
                #plt.xlim(0.9, 2.0)
                plt.title(key)
                plt.tight_layout()
            pdf.savefig(dpi=300, bbox_extra_artist=(lgd,), bbox_inches='tight')
            plt.close()
        if extra:
            figure = plt.subplot2grid((4,2), (0, 0), rowspan=2)
            imgplot = plt.imshow(img, interpolation='none')
            plt.xticks([])
            plt.yticks([])
            
            for m in range(extra):
                if chunk < 1:
                    i=0; j=0
                    key = keys[m]
                else:
                    key = keys[((plots-2)*i+j) + m+1]
                bond = bond_order_wiberg[key]
                atom_1 = molecule_charged.GetAtom(oechem.OEHasMapIdx(key[0]))
                atom_2 = molecule_charged.GetAtom(oechem.OEHasMapIdx(key[1]))
                bond_oe = molecule_charged.GetBond(atom_1, atom_2)
                wiberg_oe = bond_oe.GetData('WibergBondOrder')
                plt.subplot(plots/2, 2, 4+m)
                plt.hist(bond_order_wiberg[key], alpha=0.5, label='Wiberg_psi4')
                plt.hist(bond_order_mayer[key], alpha=0.5, label='Mayer_psi4')
                plt.vlines(wiberg_oe, ymin=0, ymax=20, label='Wiberg_OE')
                if m == 0:
                    lgd = plt.legend(prop={'size': 8}, loc='center left', bbox_to_anchor=(1, 0.5))
                #plt.xlim(0.9, 2.0)
                plt.title(key)
                plt.tight_layout()
            pdf.savefig(dpi=300, bbox_extra_artist=(lgd,), bbox_inches='tight')
            plt.close()






