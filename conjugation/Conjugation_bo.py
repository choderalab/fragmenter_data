
# coding: utf-8

# In[1]:


from IPython.display import SVG, display
from rdkit.Chem import AllChem as Chem 
from rdkit.Chem.Draw import rdMolDraw2D, MolToFile #IPythonConsole #Needed to show molecules
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions
import pandas as pd
from cairosvg import svg2png
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.image as mpimg
import seaborn as sbn

from openeye import oechem, oequacpac, oeomega, oedepict
from fragmenter import utils, fragment
import oenotebook as oenb

from openmoltools import openeye
import os
import glob
import json
from cclib.parser import Psi
from cclib.parser.utils import convertor


# In[2]:


# Generate OpenEye molecule image with Wiberg Bond order and conjugated bonds highlighted. 
# Load molecules from smi file
mollist = utils.to_oemol('kinase_inhibitors.smi', title='title', verbose=False)


# In[6]:


for molecule in mollist:
    # Generate tautomers for OpenEye conjugation model
    tautomers = fragment.expand_states(molecule, protonation=False, tautomers=True, stereoisomers=False, 
                                   suppress_hydrogens=True, carbon_hybridization=True, max_states=100, return_molecules=True)
    # Generate figure
    atom_indices = utils.tag_conjugated_bond(molecule, tautomers=tautomers)
    utils.depict_conjugation(molecule, height=700, width=1000, fname='images/{}_oe_conj.png'.format(molecule.GetTitle()), label=None)
    


# In[9]:


# Add OpenEye WBO to depiction
for molecule in mollist:
    # Generate charges
    charged = openeye.get_charges(molecule)
    charged.SetTitle(molecule.GetName())
    atom_indices = utils.tag_conjugated_bond(charged, tag='WibergBondOrder', threshold=1.05)
    utils.depict_conjugation(charged, height=700, width=1000, 
                             fname='images/{}_oe_labeled_1.05.png'.format(molecule.GetTitle()), label='WibergBondOrder')
    atom_indices = utils.tag_conjugated_bond(charged, tag='WibergBondOrder', threshold=1.2)
    utils.depict_conjugation(charged, height=700, width=1000, 
                             fname='images/{}_oe_labeled_1.2.png'.format(molecule.GetTitle()), label='WibergBondOrder')
    


# In[8]:


for molecule in mollist:
    name = molecule.GetTitle()
    print(name)
    output_files = glob.glob('geometry_opt/{}/*.output.bo.json'.format(name))
    if not output_files:
        print('No files found for {}'.format(name))
    bond_orders = {}
    for file in output_files:
        f = open(file, 'r')
        data = json.load(f)
        f.close()
        # Check for errors
        error = data['error']
        if error:
            print(error)
            continue
    data['bond_orders'] = utils.bond_order_from_psi4_raw_output(data['raw_output'])
    #bond_order.append(data['bond_orders'])
    bond_orders[data['return_value']] = data['bond_orders']
    
    molecule, atom_map = utils.get_atom_map(tagged_smiles=data['tagged_smiles'])

    boltzman_weighted_bo = utils.boltzman_average_bond_order(bond_orders)
    utils.bond_order_tag(molecule, atom_map, boltzman_weighted_bo)

    atom_indices = utils.tag_conjugated_bond(molecule, tag='Mayer_psi', threshold=1.05)
    # Generate image
    utils.highlight_bonds(mol_copy=molecule, label='Mayer_psi4', height=800, width=1000,
                   fname='images/{}_Mayer_1.05.png'.format(name))
    
    atom_indices = utils.tag_conjugated_bond(molecule, tag='Mayer_psi', threshold=1.2)
    # Generate image
    utils.highlight_bonds(mol_copy=molecule, label='Mayer_psi4', height=800, width=1000,
                   fname='images/{}_Mayer_1.2.png'.format(name))
    
    atom_indices = utils.tag_conjugated_bond(molecule, tag='Wiberg_psi', threshold=1.05)
    # Generate image
    utils.highlight_bonds(mol_copy=molecule, label='Wiberg_psi4', height=800, width=1000,
                   fname='images/{}_Wiberg_psi_1.05.png'.format(name))
    atom_indices = utils.tag_conjugated_bond(molecule, tag='Wiberg_psi', threshold=1.2)
    # Generate image
    utils.highlight_bonds(mol_copy=molecule, label='Wiberg_psi4', height=800, width=1000,
                   fname='images/{}_Wiberg_psi_1.2.png'.format(name))
    
    


# In[3]:


for mol in mollist:
    utils.highlight_bonds(mol, conjugation=False, rotor=True, width=1000, height=800, 
                          fname='images/{}_oe_rotor.png'.format(mol.GetTitle()))


# In[5]:


with PdfPages('conjugation_models.pdf') as pdf:
    for mol in mollist:
        name = mol.GetTitle()
        print(name)
        figure = plt.subplot(2, 2, 1)
        img = mpimg.imread('images/{}_rdkit_res.png'.format(name))
        imgplot = plt.imshow(img, interpolation='none')
        plt.xticks([])
        plt.yticks([])
        plt.title('RDKit Conjugation')
        
        figure = plt.subplot(2, 2, 2)
        img = mpimg.imread('images/{}_oe_conj.png'.format(name))
        imgplot = plt.imshow(img, interpolation='none')
        plt.xticks([])
        plt.yticks([])
        plt.title('OpenEye from tautomers')
        
        igure = plt.subplot(2, 2, 3)
        img = mpimg.imread('images/{}_rdkit_rotor.png'.format(name))
        imgplot = plt.imshow(img, interpolation='none')
        plt.xticks([])
        plt.yticks([])
        plt.title('RDKit rotatable bonds')
        
        figure = plt.subplot(2, 2, 4)
        img = mpimg.imread('images/{}_oe_rotor.png'.format(name))
        imgplot = plt.imshow(img, interpolation='none')
        plt.xticks([])
        plt.yticks([])
        plt.title('OpenEye rotatable bonds')
        
        plt.suptitle(name)
        plt.tight_layout()
        
        pdf.savefig(dpi=300, bbox_inches='tight')
        plt.close()
        
        figure = plt.subplot(2, 2, 1)
        img = mpimg.imread('images/{}_oe_labeled_1.2.png'.format(name))
        imgplot = figure.imshow(img, interpolation='none')
        plt.xticks([])
        plt.yticks([])
        plt.title('OE Wiberg 1.2')
        
        figure = plt.subplot(2, 2, 2)
        img = mpimg.imread('images/{}_oe_labeled_1.05.png'.format(name))
        imgplot = figure.imshow(img, interpolation='none')
        plt.xticks([])
        plt.yticks([])
        plt.title('OE Wiberg 1.05')
        
        figure = plt.subplot(2, 2, 3)
        img = mpimg.imread('images/{}_Wiberg_psi_1.05.png'.format(name))
        imgplot = figure.imshow(img, interpolation='none')
        plt.xticks([])
        plt.yticks([])
        plt.title('Psi4 Wiberg Lowdin 1.05')
        
        figure = plt.subplot(2, 2, 4)
        img = mpimg.imread('images/{}_Mayer_1.05.png'.format(name))
        imgplot = figure.imshow(img, interpolation='none')
        plt.xticks([])
        plt.yticks([])
        plt.title('Psi4 Mayer 1.05')
        
        plt.suptitle(name)
        
        plt.tight_layout()
        
        pdf.savefig(dpi=300, bbox_inches='tight')
        plt.close()
        #plt.close()


# In[20]:


import numpy as np
# First generate bond order distribution
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
        if bond.IsInRing():
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
    with PdfPages('{}_bond_orders_no_rings_2.pdf'.format(ki)) as pdf:
        n = len(bond_order_wiberg)
        plots = 6
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
                figure = plt.subplot2grid((3,2), (0, 0), rowspan=2)
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
                    plt.subplot2grid((3, 2), (0, 1))
                    #plt.(bond_order_wiberg[key], label='Wiberg_psi4', kde=True)
                    plt.hist(bond_order_wiberg[key], label='Wiberg_psi4', alpha=1.0)
                    plt.hist(bond_order_mayer[key], label='Mayer_psi4', alpha=1.0)
                    plt.vlines(wiberg_oe, ymin=0, ymax=len(bond_order_wiberg[key])/2.5, label='Wiberg_OE', linewidth=0.9)
                    plt.vlines(1.2, ymin=0, ymax=len(bond_order_wiberg[key])/2.5, label='Current threshold', 
                           linewidth=0.9, color='red')
                    #plt.xlim(0.8, 3.3)
                    #plt.ylim(0, len(bond_order_wiberg[key])/2.5)
                    plt.title(key)
                    lgd = plt.legend(prop={'size': 8}, loc='center left', bbox_to_anchor=(1, 0.5))
                #if j == 1:
                #    plt.subplot(plots/2, 2, 3)
                #if j == 2:
                #    plt.subplot(plots/2, 2, 4)
                
                else:
                    plt.subplot(plots/2, 2, j+3)
                    plt.hist(bond_order_wiberg[key], label='Wiberg_psi4', alpha=1.0)
                    plt.hist(bond_order_mayer[key], label='Mayer_psi4', alpha=1.0)
                    plt.vlines(wiberg_oe, ymin=0, ymax=len(bond_order_wiberg[key])/2.5, label='Wiberg_OE', linewidth=0.9)
                    plt.vlines(1.2, ymin=0, ymax=len(bond_order_wiberg[key])/2.5, label='Current threshold', 
                               linewidth=0.9, color='red')

                    #plt.xlim(0.8, 3.3)
                    #plt.ylim(0, len(bond_order_wiberg[key])/2.5)
                    #plt.xlim(0.9, 2.0)
                    plt.title(key)
                plt.tight_layout()
            pdf.savefig(dpi=300, bbox_extra_artist=(lgd,), bbox_inches='tight')
            plt.close()
            
        if extra:
            figure = plt.subplot2grid((3,2), (0, 0), rowspan=2)
            img = mpimg.imread('{}_mapped.png'.format(ki))
            imgplot = plt.imshow(img, interpolation='none')
            #imgplot = plt.imshow(img, interpolation='none')
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
                if m == 0:
                    plt.subplot(plots/2, 2, 2)
                else:
                    plt.subplot(plots/2, 2, m+3)
                #plt.subplot(plots/2, 2, 3+m)
                plt.hist(bond_order_wiberg[key], label='Wiberg_psi4', alpha=1.0)
                plt.hist(bond_order_mayer[key],  label='Mayer_psi4', alpha=1.0)
                plt.vlines(wiberg_oe, ymin=0, ymax=len(bond_order_wiberg[key])/2.5, label='Wiberg_OE', linewidth=0.9)
                plt.vlines(1.2, ymin=0, ymax=len(bond_order_wiberg[key])/2.5, label='Current threshold', 
                           linewidth=0.9, color='red')
                #plt.xlim(0.8, 3.3)
                #plt.ylim(0, len(bond_order_wiberg[key])/2.5)
                if m == 0:
                    lgd = plt.legend(prop={'size': 8}, loc='center left', bbox_to_anchor=(1, 0.5))
                #plt.xlim(0.9, 2.0)
                plt.title(key)
                plt.tight_layout()
            pdf.savefig(dpi=300, bbox_extra_artist=(lgd,), bbox_inches='tight')
            plt.close()

