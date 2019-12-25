import json
import numpy as np
import os
import warnings

from openeye import oechem
from fragmenter import chemi, utils

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib.colors as cl
import matplotlib as mpl
from matplotlib import gridspec
import seaborn as sbn

def mmd_x_xsqred(x, y):
    """
    Maximum mean discrepancy with squared kernel
    This will distinguish mean and variance
    see https://stats.stackexchange.com/questions/276497/maximum-mean-discrepancy-distance-distribution
    Parameters
    ----------
    x : list of ints
    y : list of ints

    Returns
    -------
    mmd score

    """

    y_arr = np.asarray(y)
    y_squared = y_arr*y_arr
    x_arr = np.asarray(x)
    x_squared = np.square(x_arr)

    E_x = np.mean(x_arr)
    E_y = np.mean(y_arr)

    E_x_squared = np.mean(x_squared)
    E_y_squared = np.mean(y_squared)

    mmd2 = (E_x - E_y)**2 + (E_x_squared - E_y_squared)**2
    return mmd2

def get_bond(mol, bond_idx):
    a1 = mol.GetAtom(oechem.OEHasMapIdx(bond_idx[0]))
    a2 = mol.GetAtom(oechem.OEHasMapIdx(bond_idx[1]))
    bond = mol.GetBond(a1, a2)
    if not bond:
        raise ValueError("({}) atoms are not connected".format(bond_idx))
    return bond

def find_parent_smiles(wbo_dict, mapped=True):
    for bond in wbo_dict:
        for frag in wbo_dict[bond]:
            if 'parent' in frag:
                if mapped:
                    return wbo_dict[bond][frag]['map_to_parent']
                else:
                    return frag
                #mapped_smiles = wbo_dict[bond][frag]['map_to_parent']
                #return mapped_smiles


def sort_fragments(frags):
    sorted_frags = sorted([(frags[f]['elf_estimate'], f) for f in frags])
    return [f[1] for f in sorted_frags]

def find_highest_score(full_frags):
    s = 0
    for b in full_frags:
        for f in full_frags[b]:
            if full_frags[b][f]['mmd_exp'] > s:
                s = full_frags[b][f]['mmd_exp']
    return s

def rbg_to_int(rbg, alpha):
    rbg[-1] = int(rbg[-1]*alpha)
    colors = [int(round(i*255)) for i in rbg[:-1]]
    colors.append(int(rbg[-1]))
    return colors

def generate_molecule_images(fragments_dict, name, colors, parent_mol):
    for i, b in enumerate(fragments_dict):
        to_plot = []
        wbos = []
        sorted_frags = sort_fragments(fragments_dict[b])
        for f in sorted_frags:
            smiles = fragments_dict[b][f]['map_to_parent']
            mol = oechem.OEMol()
            oechem.OESmilesToMol(mol, smiles)
            mol.SetTitle(str(fragments_dict[b][f]['mmd_exp']))
            to_plot.append(mol)
            wbos.append(fragments_dict[b][f]['elf_estimate'])
        int_colors = [rbg_to_int(rbg, alpha=150) for rbg in colors[i]]
        colors_oe = [oechem.OEColor(*j) for j in int_colors]
        fname = 'selected/{}/rescore/{}_bond_{}_{}_aligned.pdf'.format(name, name, str(b[0]), str(b[1]))
        chemi.to_pdf(to_plot, fname, rows=3, cols=3, bond_map_idx=b, bo=wbos, color=colors_oe,
                     align=to_plot[0])
        fname = 'selected/{}/rescore/{}_bond_{}_{}_aligned_to_parent.pdf'.format(name, name, str(b[0]), str(b[1]))
        chemi.to_pdf(to_plot, fname, rows=3, cols=3, bond_map_idx=b, bo=wbos, color=colors_oe,
                     align=parent_mol)
        fname = 'selected/{}/rescore/{}_bond_{}_{}_no_alignemnet.pdf'.format(name, name, str(b[0]), str(b[1]))
        chemi.to_pdf(to_plot, fname, rows=3, cols=3, bond_map_idx=b, bo=wbos, color=colors_oe)

def fragment_wbo_ridge_plot(data, filename, rug=True):
    """
    data: dict of dict
        bond: frags: [wbos]
    colors: numpy array
        output from matplotlib color map
    sorted_keys: dict
        bond: sorted_frags
    filname: str
        filename
    """
    higlight_for_bonds = []
    with PdfPages(filename) as pdf:
        for bond in data:
            # sort fragments by estimated wbo
            sorted_frags = sort_fragments(data[bond])
            n = len(data[bond])
            fig = plt.figure()
            gs = gridspec.GridSpec(n, 2, width_ratios=[20, 1])
            fig.dpi = 400
            x_min = 3
            x_max = 0
            for f in data[bond]:
                wbo = data[bond][f]['individual_confs']
                if min(wbo) < x_min:
                    x_min = min(wbo)
                if max(wbo) > x_max:
                    x_max = max(wbo)

            scores = [data[bond][f]['mmd_exp'] for f in sorted_frags]
            norm = plt.Normalize(0, max(scores))
            colors = plt.cm.plasma(norm(scores))
            higlight_for_bonds.append(colors)
            for i, frag in enumerate(sorted_frags):
                wbo = data[bond][frag]['individual_confs']

                wbo_s = data[bond][frag]['elf_estimate']
                ax = ax1 = plt.subplot(gs[i, :-1])
                ax.spines['left'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.spines['top'].set_visible(False)
                ax.spines['bottom'].set_visible(False)
                ax.patch.set_facecolor('none')
                if 'parent' in frag:
                    sbn.kdeplot(wbo, shade= True, color='cyan', alpha=1.0)
                    sbn.kdeplot(wbo, lw=1.5, color='black')
                else:
                    sbn.kdeplot(wbo, shade=True, color=colors[i], alpha=0.3)
                    sbn.kdeplot(wbo, lw=0.4, color=colors[i])
                if rug:
                    sbn.distplot(wbo, hist=False, rug=True, kde=False, color='black')

                # img = plt.axvline(x=wbo_s, ymin=0, ymax=1, color='black', linewidth=0.5)
                if len(wbo) < 2:
                    plt.axvline(x=wbo_s, ymin=0, ymax=1, color=colors[i], linewidth=2.0, alpha=0.8)
                if rug:
                    plt.axvline(x=wbo_s, ymin=0, ymax=1, color='black', linewidth=0.5)
                plt.xlim(x_min-0.1, x_max+0.1)
                plt.yticks([])
                ax.yaxis.set_label_coords(-0.05, 0)
                if i != n-1:
                    plt.xticks([])
                else:
                    plt.xlabel('Bond order')
                if i == 0:
                    plt.title(bond)

            # Add colorbar
            ax = plt.subplot(gs[:, -1])
            cmap = mpl.cm.plasma
            norm = mpl.colors.Normalize(min(scores), max(scores))
            cb1 = mpl.colorbar.ColorbarBase(ax, cmap=cmap,
                                norm=norm,
                                orientation='vertical')
            cb1.set_label('mmd_score')

            # Magic to get overlapping distributions
            overlap=0.5
            h_pad = 5 + (- 5*(1 + overlap))
            fig.tight_layout(h_pad=h_pad)
            pdf.savefig(bbox_inches='tight')
            plt.close()
        return higlight_for_bonds


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="calculate OE WBO")
    parser.add_argument('-n', '--name', type=str, help='Molecule name with number of its state')
    args = parser.parse_args()
    name = args.name

    infile = '../fragment_bond_orders/validation_set/{}/{}_oe_wbo_by_bond.json'.format(name, name)
    try:
        with open(infile, 'r') as f:
            bonds_dist = json.load(f)
    except:
        filename = '{}_one_parent_map.json'.format(infile.split('.json')[0])
        with open(filename, 'r') as f:
            bonds_dist = json.load(f)
    bonds_file = 'selected/{}/{}_selected_bonds.json'.format(name, name)
    with open(bonds_file, 'r') as f:
        bonds = json.load(f)['bonds']

    try:
        os.mkdir('selected/{}/rescore/'.format(name))
    except FileExistsError:
         print('{} directory already exists. Files will be overwritten'.format(name))

    parent_mapped_smiles = find_parent_smiles(bonds_dist, mapped=True)
    # Generate image with map labels
    parent_mol = oechem.OEMol()
    oechem.OESmilesToMol(parent_mol, parent_mapped_smiles)
    chemi.mol_to_image_atoms_label(parent_mol, fname='selected/{}/{}.png'.format(name, name))


    # Only keep fragments if all 1-5 atoms are around central bonds
    full_frags = {}
    for bond in bonds:
        b = get_bond(parent_mol, bond)
        if not b:
            warnings.warn('bond {} does not exist in {}'.format(bond, name))

        bond_key = tuple(bond)
        serialized_key = utils.serialize_bond(bond_key)
        full_frags[bond_key] = bonds_dist[serialized_key]

    # Calculate mmd score for each fragment
    parent_key = find_parent_smiles(full_frags, mapped=False)
    for b in full_frags:
        parent_dist = full_frags[b][parent_key]['individual_confs']
        for f in full_frags[b]:
            mmd = mmd_x_xsqred(x=parent_dist, y=full_frags[b][f]['individual_confs'])
            full_frags[b][f]['mmd_exp'] = mmd

    # serialize and save
    serialized_full_frags = {}
    for b in full_frags:
        b_key = utils.serialize_bond(b)
        serialized_full_frags[b_key] = full_frags[b]

    with open('selected/{}/rescore/{}_oe_wbo_with_score.json'.format(name, name), 'w') as f:
            json.dump(serialized_full_frags, f, indent=2, sort_keys=True)

    # Save scores separately - this will make it easier to rank and finalize validation set
    # Save all scores
    scores = {}
    frags_with_scores = {}
    for b in full_frags:
        b_key = utils.serialize_bond(b)
        scores[b_key] = []
        frags_with_scores[b_key] = {'frags': [], 'mmd_scores': []}
        for f in full_frags[b]:
            scores[b_key].append(full_frags[b][f]['mmd_exp'])
            frags_with_scores[b_key]['frags'].append(f.split('_')[0])
            frags_with_scores[b_key]['mmd_scores'].append(full_frags[b][f]['mmd_exp'])

    scores['greatest_discrepancy'] = find_highest_score(full_frags)

    with open('selected/{}/rescore/{}_mmd_exp_scores.json'.format(name, name), 'w') as f:
        json.dump(scores, f, indent=2, sort_keys=True)
    with open('selected/{}/rescore/{}_frag_with_scores.json'.format(name, name), 'w') as f:
        json.dump(frags_with_scores, f, indent=2, sort_keys=True)

    # Plot joyplot of fragment distribution sorted by elf wbo and colored by mmd score
    colors = fragment_wbo_ridge_plot(full_frags, filename='selected/{}/rescore/{}_ridge_with_rug.pdf'.format(name, name))
    colors = fragment_wbo_ridge_plot(full_frags, filename='selected/{}/rescore/{}_ridge.pdf'.format(name, name), rug=False)

    # Generate images of molecules
    generate_molecule_images(full_frags, name, colors, parent_mol)
