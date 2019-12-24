"""
This script generates figures to illustrate what the results of the combinatorial fragmentation looks like.
I chose Debrafenib bond (35, 14) because it illustrates really well how specific important chemical cuts can shift the
WBO distributions even if the fragment has all other parts of the parent molecules. Most other results are a bit
more complicated to explain as clearly or they don't have such nice 3/4 distinct clusters.
The first figure shows all distributions of WBOs for all fragments that contain the bond (35, 14) of Dabrafenib
The second figure is the Perato front
The third figure generates selected fragments to illustrate the points mentioned above
"""

import fragmenter
from openeye import oechem, oedepict, oegraphsim

import json
import numpy as np
import matplotlib as mpl
import seaborn as sbn
import matplotlib.pyplot as plt


def sort_fragments_by_elf10wbo(frags):
    """
    Sort fragments by ELF10 WBO. This helps with plotting all distributions so the distributions of the different
    clusters are together
    Parameters
    ----------
    frags : dict
        {'smiles': {'ensamble': , 'individual_confs': []}

    Returns
    -------
    list of sorted frags in tuple sorted by elf10 wbos [(elf10 wbo, smiles)]
    """
    sorted_frags = sorted([(frags[f]['ensamble'], f) for f in frags])
    return sorted_frags


def rbg_to_int(rbg, alpha):
    """
    Convert rbg color to ints for openeye
    Parameters
    ----------
    rbg : list
        rbg
    alpha : int

    Returns
    -------
    list of ints

    """
    rbg[-1] = int(rbg[-1]*alpha)
    colors = [int(round(i*255)) for i in rbg[:-1]]
    colors.append(int(rbg[-1]))
    return colors


def n_heavy_atoms(smiles):
    """
    Count number of heavy atoms in SMILES
    Parameters
    ----------
    smiles: str
        SMILES to count heavy atoms
    Returns
    -------
    n: int
        number of heavy atoms in molecules
    """
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, smiles)
    n = 0
    for a in mol.GetAtoms():
        if not a.IsHydrogen():
            n += 1
    return n


def get_bond(mol, bond_idx):
    """
    Get oebond object
    Parameters
    ----------
    mol : oemol
        Molecule to extract bond from
    bond_idx : tuple of ints
        tuple of map indices of atoms in bond

    Returns
    -------
    bond: oebond

    """
    atoms = [mol.GetAtom(oechem.OEHasMapIdx(i)) for i in bond_idx]
    bond = mol.GetBond(atoms[0], atoms[1])
    return bond


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
    return np.sqrt(mmd2)


def fragment_wbo_ridge_plot(data, filename, indices=None):
    """
    Generate ridge plot of all WBO distribtuions
    data: dict
        {frag: {'mmd_score': float, 'ensamble': float, 'individual_conf': [wbos]}
    filname: str
        filename
    indices: list of ints.
        This is used to find the distributions to outline
    Returns
    -------
    highlights: np.array
        results from matplotlib normed colormap.
        RBG listed in order of sorted frags
    """
    higlight_for_bonds = []

    # sort fragments by estimated wbo
    sorted_frags = sort_fragments_by_elf10wbo(data)
    sorted_frags = [f[1] for f in sorted_frags]
    n = len(data)
    fig = plt.figure()
    gs = mpl.gridspec.GridSpec(n, 2, width_ratios=[20, 1])
    fig.dpi = 400
    x_min = 3
    x_max = 0
    for f in data:
        wbo = data[f]['individual_conf']
        if min(wbo) < x_min:
            x_min = min(wbo)
        if max(wbo) > x_max:
            x_max = max(wbo)

    scores = [data[f]['mmd_score'] for f in sorted_frags]
    norm = plt.Normalize(0, max(scores))
    colors = plt.cm.viridis_r(norm(scores))
    higlight_for_bonds.append(colors)
    for i, frag in enumerate(sorted_frags):
        wbo = data[frag]['individual_conf']

        wbo_s = data[frag]['ensamble']
        ax = ax1 = plt.subplot(gs[i, :-1])
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.patch.set_facecolor('none')
        if 'parent' in frag:
            sbn.kdeplot(wbo, shade= True, color='red', alpha=1.0)
            sbn.kdeplot(wbo, lw=2, color='black')

        else:
            sbn.kdeplot(wbo, shade=True, color=colors[i], alpha=0.2)
            sbn.kdeplot(wbo, lw=0.3, color=colors[i])
            if i in indices:
                sbn.kdeplot(wbo, shade=True, color=colors[i], alpha=0.3)
                sbn.kdeplot(wbo, lw=1.5, color='black')

        if len(wbo) < 2:
            xs = np.asarray([wbo[0]-0.0001, wbo[0], wbo[0]+0.0001])
            sbn.kdeplot(xs, shade=True, color=colors[i], alpha=0.3)
            sbn.kdeplot(xs, lw=0.3, color=colors[i])
            if i in indices:
                 sbn.kdeplot(xs, lw=1, color='black')

        plt.xlim(x_min-0.02, x_max+0.02)
        plt.yticks([])
        ax.yaxis.set_label_coords(-0.05, 0)
        if i != n-1:
            plt.xticks([])
        else:
            plt.xlabel('Bond order')
        #if i == 0:
        #    plt.title(bond)

    # Add colorbar
    ax = plt.subplot(gs[:, -1])
    cmap = mpl.cm.viridis_r
    norm = mpl.colors.Normalize(min(scores), max(scores))
    cb1 = mpl.colorbar.ColorbarBase(ax, cmap=cmap,
                        norm=norm,
                        orientation='vertical')
    cb1.set_label('Distance score')

    # Magic to get overlapping distributions
    overlap=0.5
    h_pad = 5 + (- 5*(1 + overlap))
    fig.tight_layout(h_pad=h_pad)
    plt.savefig(fname=filename, bbox_inches='tight')
    plt.close()
    return higlight_for_bonds


def plot_pareto_frontier(Xs, Ys, parent_idx, maxX=True, maxY=True, filename=None, title=None, indices=None):
    """
    Plot Pareto front

    Parameters:
    ----------
    Xs: list of floats
    Ys: list of floats
    maxX: bool
        Used to know where to plot the front. If True, at maxX, if False, at minX
    xaxY: bool
        Used to know where to plot front. If True, maxY, if False, minY
    filename: str
    title: str
    indices: list of ints, optional, default None
        If given, these data points will be highlighted and annotated

    Returns:
    -------
    front: list
        List of (score, size) of data points on front

    """
    '''Pareto frontier selection process'''
    sorted_list = sorted([[Xs[i], Ys[i]] for i in range(len(Xs))], reverse=maxY)
    pareto_front = [sorted_list[0]]
    for i, pair in enumerate(sorted_list[1:]):
        if maxY:
            if pair[1] >= pareto_front[-1][1]:
                pareto_front.append(pair)

        else:
            if pair[1] <= pareto_front[-1][1]:
                pareto_front.append(pair)

    '''Plotting process'''
    fig, ax = plt.subplots()
    # ax.scatter(Xs,Ys,edgecolors='black', c='white')

    pf_X = [pair[0] for pair in pareto_front]
    pf_Y = [pair[1] for pair in pareto_front]
    plt.plot(pf_X, pf_Y, color='black', linewidth=0.8, zorder=1)

    for i, j in enumerate(indices):
        if j == parent_idx:
            ax.scatter(Xs[j], Ys[j], marker='o', edgecolors='red', c='red', s=60, zorder=2)
            ax.annotate(i + 1, (Xs[j], Ys[j]), textcoords='offset points', xytext=(5, 4), color='red', fontsize=7)
            ax.annotate('   (parent)', (Xs[j], Ys[j]), textcoords='offset points', xytext=(5, 4), color='red',
                        fontsize=7)
        else:
            ax.scatter(Xs[j], Ys[j], marker='o', edgecolors='red', c='black', s=60, zorder=2)
            ax.annotate(i + 1, (Xs[j], Ys[j]), textcoords='offset points', xytext=(5, 4), color='black', fontsize=7)
    img = ax.scatter(Xs, Ys, c=Xs, edgecolors='black', cmap=plt.cm.get_cmap('viridis_r'), alpha=1,
                     zorder=3)

    fig.colorbar(img, ax=ax)

    plt.plot(pf_X, pf_Y, '.', color='black', markersize=2, zorder=4)

    plt.xlabel("MMD score")
    plt.ylabel("Computational cost (molecular size)^3")
    plt.xlim(-0.005, max(Xs) + 0.01)
    plt.title(title)
    if filename:
        plt.savefig(filename)
    plt.close()
    return pareto_front


def visualize_mols(smiles, fname, rows, cols, bond_idx, wbos, colors, align_to=0):
    """
    Visualize molecules with highlighted bond and labeled with WBO
    Parameters
    ----------
    smiles : list of SMILES to visualize.
        bond atoms should have map indices
    fname : str
        filename
    rows : int
    cols : int
    bond_idx : tuple of atom maps of bond to highlight.
    wbos : list of floats
    colors : list of hex values for colors
    align_to: int, optional, default 0
        index for which molecule to align to. If zero, will align to first molecules in SMILES list

    """
    itf = oechem.OEInterface()

    ropts = oedepict.OEReportOptions(rows, cols)
    ropts.SetHeaderHeight(25)
    ropts.SetFooterHeight(25)
    ropts.SetCellGap(2)
    ropts.SetPageMargins(10)
    report = oedepict.OEReport(ropts)

    cellwidth, cellheight = report.GetCellWidth(), report.GetCellHeight()
    opts = oedepict.OE2DMolDisplayOptions(cellwidth, cellheight, oedepict.OEScale_AutoScale)
    oedepict.OESetup2DMolDisplayOptions(opts, itf)

    # align to chosen molecule
    ref_mol = oechem.OEGraphMol()
    oechem.OESmilesToMol(ref_mol, smiles[align_to])
    oedepict.OEPrepareDepiction(ref_mol)

    mols = []
    for s in smiles:
        mol = oechem.OEMol()
        oechem.OESmilesToMol(mol, s)
        mols.append(mol)
        oedepict.OEPrepareDepiction(mol, False, True)
        minscale = float("inf")
        minscale = min(minscale, oedepict.OEGetMoleculeScale(mol, opts))

    opts.SetScale(minscale)
    for i, mol in enumerate(mols):

        cell = report.NewCell()
        oedepict.OEPrepareDepiction(mol, False, True)
        bond = get_bond(mol, bond_idx)
        atom_bond_set = oechem.OEAtomBondSet()
        atom_bond_set.AddAtoms([bond.GetBgn(), bond.GetEnd()])
        atom_bond_set.AddBond(bond)

        hstyle = oedepict.OEHighlightStyle_BallAndStick
        if i == 3:
            hcolor = oechem.OERed
        else:
            hcolor = oechem.OEColor(colors[i])

        overlaps = oegraphsim.OEGetFPOverlap(ref_mol, mol, oegraphsim.OEGetFPType(oegraphsim.OEFPType_Tree))
        oedepict.OEPrepareMultiAlignedDepiction(mol, ref_mol, overlaps)

        disp = oedepict.OE2DMolDisplay(mol, opts)
        oedepict.OEAddHighlighting(disp, hcolor, hstyle, atom_bond_set)

        bond_label = oedepict.OEHighlightLabel("{:.2f}".format((wbos[i])), hcolor)
        oedepict.OEAddLabel(disp, bond_label, atom_bond_set)
        oedepict.OERenderMolecule(cell, disp)
        # oedepict.OEDrawCurvedBorder(cell, oedepict.OELightGreyPen, 10.0)

    return (oedepict.OEWriteReport(fname, report))


if __name__ == '__main__':
    with open('Dabrafenib_oe_wbo_by_bond_with_parent.json', 'r') as f:
        wbos = json.load(f)
    with open('Dabrafenib_fragments.json', 'r') as f:
        fragments = json.load(f)

    # Deserialize
    wbos_des = {}
    for bond in wbos:
        bond_des = fragmenter.utils.deserialize_bond(bond)
        wbos_des[bond_des] = wbos[bond]

    # Calculate scores and store
    for bond in wbos_des:
        parent = wbos_des[bond]['parent']
        for frag in wbos_des[bond]:
            mmd_score = mmd_x_xsqred(wbos_des[bond][frag]['individual_conf'], parent['individual_conf'])
            wbos_des[bond][frag]['mmd_score'] = mmd_score

    # Indices of fragments to illustrate in figure
    indices = [0, 1, 2, 7, 10, 12, 18, 20, 26, 34, 36]
    highlights = fragment_wbo_ridge_plot(wbos_des[(35, 14)], 'ridge_plot.pdf', indices=indices)

    # Prepare data for Perato front.
    sorted_elf10wbo_frags = sort_fragments_by_elf10wbo(wbos_des[(35, 14)])
    sorted_frags = [i[1] for i in sorted_elf10wbo_frags]
    sorted_scores = [wbos_des[(35, 14)][f]['mmd_score'] for f in sorted_frags]

    parent_smiles = "[H:36][c:1]1[c:3]([c:9]([c:13]([c:10]([c:4]1[H:39])[N:28]([H:55])[S:35](=[O:29])(=[O:30])[c:14]2[c:11]([c:5]([c:2]([c:6]([c:12]2[F:32])[H:41])[H:37])[H:40])[F:31])[F:33])[c:16]3[c:17]([s:34][c:18]([n:25]3)[C:23]([C:20]([H:44])([H:45])[H:46])([C:21]([H:47])([H:48])[H:49])[C:22]([H:50])([H:51])[H:52])[c:15]4[c:7]([c:8]([n:24][c:19]([n:26]4)[N:27]([H:53])[H:54])[H:43])[H:42])[H:38]"

    sorted_map_to_parent = []
    for frag in sorted_frags:
        if frag == 'parent':
            sorted_map_to_parent.append(parent_smiles)
            continue
        map_to_parent = fragments[frag]['provenance']['routine']['enumerate_fragments']['map_to_parent']
        sorted_map_to_parent.append(map_to_parent)
    frag_sizes = [n_heavy_atoms(s)**3 for s in sorted_map_to_parent]

    front = plot_pareto_frontier(sorted_scores, frag_sizes, maxX=False, maxY=False, indices=indices, parent_idx=7,
                                 filename='pareto_front.pdf')

    # Prepare data for molecules to visualize
    to_plot = []
    wbos = []
    smiles_to_plot = []
    for i in indices:
        f = sorted_frags[i]
        map_to_parent = sorted_map_to_parent[i]
        smiles_to_plot.append(map_to_parent)
        wbos.append(wbos_des[(35, 14)][f]['ensamble'])
    int_colors = [rbg_to_int(rbg, alpha=250) for rbg in highlights[0]]
    colors_oe = [oechem.OEColor(*int_colors[j]) for j in indices]

    visualize_mols(smiles_to_plot, 'fragments.pdf', rows=3, cols=2, wbos=wbos, colors=colors_oe, bond_idx=(35, 14),
                   align_to=6)

