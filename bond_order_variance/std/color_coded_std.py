
# coding: utf-8

# In[1]:


#get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
from fragmenter import chemi
from openeye import oechem
import oenotebook as oenb
import json
import numpy as np
import glob
import cmiles
import os
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sbn
import pandas as pd
import matplotlib.image as mpimg


# In[19]:


oemols = chemi.file_to_oemols('../../combinatorial_fragmentation/kinase_inhibitors.smi')


# In[168]:

for mol in oemols:
    name = mol.GetTitle()
    print(name)
    output = []
    output_files = glob.glob('../../conjugation/geometry_opt/{}/*.bo.json'.format(name))
    if not output_files:
        continue
    for file in output_files:
        with open(file, 'r') as f:
            data = json.load(f)
        error = data['error']
        if error:
            print(error)
            continue
        data['bond_orders'] = chemi.bond_order_from_psi4_raw_output(data['raw_output'])
        data.pop('raw_output')
        output.append(data)


    # In[169]:


    def collect_all_std(data, hydrogen=False, rings=False, only_rings=False, halogens=False, carbonyls=False, nitriles=False):

        bond_order_std = {}
        conformations = len(data)
        mol = cmiles.utils.load_molecule(data[0]['tagged_smiles'])
        n_atoms = mol.NumAtoms()
        bond_order_wiberg = {}
        bond_order_mayer = {}
        for bond in mol.GetBonds():
            if only_rings and not rings:
                raise RuntimeError("If only rings is True, rings must be true")
            if only_rings:
                if not bond.IsInRing():
                    continue
            if not rings:
                if bond.IsInRing():
                    continue
            atom_1 = bond.GetBgn()
            atom_2 = bond.GetEnd()
            if not hydrogen:
                if atom_1.IsHydrogen() or atom_2.IsHydrogen():
                    continue
            if not halogens:
                if atom_1.IsHalogen() or atom_2.IsHalogen():
                    continue
            if not nitriles:
                if bond.GetOrder() == 3:
                    continue
            if not carbonyls:
                if bond.GetOrder() == 2 and (atom_1.IsOxygen() or atom_2.IsOxygen()):
                    continue
            map_1 = atom_1.GetMapIdx()
            map_2 = atom_2.GetMapIdx()
            bond_order_wiberg[(map_1, map_2)] = np.zeros(conformations)
            bond_order_mayer[(map_1, map_2)] = np.zeros(conformations)

        # Populate array
        for k, d in enumerate(data):
            wiberg = d['bond_orders']['Wiberg_psi4']
            mayer = d['bond_orders']['Mayer_psi4']
            for i, j in bond_order_wiberg:
                bond_order_wiberg[(i, j)][k] = wiberg[i-1][j-1]
                bond_order_mayer[(i, j)][k] = mayer[i-1][j-1]
        bond_order_std['wiberg_bo'] = bond_order_wiberg
        bond_order_std['mayer_bo'] = bond_order_mayer
        # calculate variance
        bonds = []
        mayer_std = []
        wiberg_std = []
        for bond in bond_order_mayer:
            bonds.append(bond)
            mayer_std.append(np.std(bond_order_mayer[bond]))
            wiberg_std.append(np.std(bond_order_wiberg[bond]))
        bond_order_std['wiberg_std'] = wiberg_std
        bond_order_std['mayer_std'] = mayer_std
        bond_order_std['bonds'] = bonds
        return bond_order_std
    bond_order_std_of_interest = collect_all_std(output)
    bond_order_std_rings = collect_all_std(output, only_rings=True, rings=True)
    # all bond orders
    bond_order_std = collect_all_std(output, rings=True, hydrogen=True, halogens=True, nitriles=True, carbonyls=True)
    bond_order_std_no_rings = collect_all_std(output, rings=False, hydrogen=False, halogens=True, nitriles=True, carbonyls=True)
    bond_order_std_no_h = collect_all_std(output, rings=True, hydrogen=False, halogens=True, nitriles=True, carbonyls=True)


    # In[170]:


    bond_order_std_of_interest['wiberg_bo'].keys()


    # In[171]:


    len(bond_order_std_of_interest['wiberg_bo'])


    # In[172]:


    bonds = list(bond_order_std_of_interest['wiberg_bo'].keys())
    colors = chemi._KELLYS_COLORS
    n = len(bond_order_std_of_interest['wiberg_bo'])
    fig, axes = plt.subplots(n, 1)
    fig.dpi = 400
    x_min = 3
    x_max = 0
    for b in bond_order_std_of_interest['wiberg_bo']:
        wbo = bond_order_std_of_interest['wiberg_bo'][b]
        if min(wbo) < x_min:
            x_min = min(wbo)
        if max(wbo) > x_max:
            x_max = max(wbo)

    for i, bond in enumerate(bonds):
        wbo = bond_order_std_of_interest['wiberg_bo'][bond]
        ax = plt.subplot(n, 1, i+1)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.patch.set_facecolor('none')
        sbn.kdeplot(wbo, shade= True, alpha=0.85, color=colors[i])
        #sbn.distplot(wbo, hist=False, rug=True, kde=False, color='black')
        sbn.kdeplot(wbo, lw=1.5, color=colors[i])
        #plt.axvline(x=wbo_s, ymin=0, ymax=1, color='black', linewidth=0.5)
        plt.xlim(x_min-0.05, x_max+0.05)
        plt.yticks([])
        #ax.set_yticklabels(bond_order_std_rot_bonds_opt['Imatinib']['bonds'], fontsize=6, rotation=0)
        #ax.yaxis.grid(False)
        ax.yaxis.set_label_coords(-0.05, 0)
        plt.ylabel(bond, rotation=0, size=8)
        if i != n-1:
            plt.xticks([])
        else:
            plt.xlabel('Bond order')
        if i == 0:
            #plt.legend(prop={'size': 10}, bbox_to_anchor=(1.35, 1))
            plt.title('{} WBO'.format(name))
    overlap=0.5
    h_pad = 5 + (- 5*(1 + overlap))
    fig.tight_layout(h_pad=h_pad)
    plt.savefig('{}_single_bond_wbo.pdf'.format(name))


    # In[173]:


    chemi.highlight_bond_by_map_idx(data['tagged_smiles'], bonds, fname='{}_highlighted_bonds.pdf'.format(name),
                                    map_idx=True, label_scale=1)

