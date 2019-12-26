"""
This script finds if the fragment resulting from the benchmark schemes fall on the Pareto front
"""

import fragmenter
from openeye import oechem, oedepict, oegraphsim
import json

import matplotlib.pyplot as plt
import numpy as np


def rbg_to_int(rbg, alpha):
    rbg[-1] = int(rbg[-1]*alpha)
    colors = [int(round(i*255)) for i in rbg[:-1]]
    colors.append(int(rbg[-1]))
    return colors

def n_heavy_atoms(smiles):
    """
    """
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, smiles)
    n = 0
    for a in mol.GetAtoms():
        if not a.IsHydrogen():
            n += 1
    return n

def get_bond(mol, bond_tuple):
    """
    Get bond in molecule
    Parameters
    ----------
    mol : oemole with map indices
    bond_tuple : tuple with map indices of bond

    Returns
    -------
    bond if found, False otherwise
    """

    a1 = mol.GetAtom(oechem.OEHasMapIdx(bond_tuple[0]))
    a2 = mol.GetAtom(oechem.OEHasMapIdx(bond_tuple[1]))
    if not a1 or not a2:
        return False
    bond = mol.GetBond(a1, a2)
    if not bond:
        return False
    return bond

def get_bond_nbr(mol, bond):
    # Find all 1-4 atoms around central bond
    nbrs = set()
    a1 = bond.GetBgn()
    a2 = bond.GetEnd()
    for a in a1.GetAtoms():
        if a.IsHydrogen():
            continue
        nbrs.add(a)
    for a in a2.GetAtoms():
        if a.IsHydrogen():
            continue
        nbrs.add(a)
    #nbrs = set()
  #  for nbr in nbrs_1:
  #      # keep all nbrs of nbrs
  #      for a in nbr.GetAtoms():
  #          if a.IsHydrogen():
  #              continue
  #          nbrs.add(a)
    # convert to map idxs
    nbrs_map_idx = set()
    for a in nbrs:
        nbrs_map_idx.add(a.GetMapIdx())

    return list(nbrs_map_idx)

def frags_with_nbrs(frags, nbrs):
    return_frags = {}
    for f in frags:
        keep = True
        mol = oechem.OEMol()
        oechem.OESmilesToMol(mol, frags[f]['map_to_parent'])
        for n in nbrs:
            a = mol.GetAtom(oechem.OEHasMapIdx(n))
            if not a:
                keep = False
        if keep:
            return_frags[f] = frags[f]
    return return_frags

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


def to_pdf(molecules, bond_map_idx, fname, rows=3, cols=2, align=None, colors=None, wbos=None, names=None):
    """
    Generate PDF of list of oemols or SMILES

    Parameters
    ----------
    molecules : list of OEMols or SMILES
    fname : str
        Name of PDF
    rows : int
        How many rows of molecules per page
    cols : int
        How many columns of molecule per page
    bond_map_idx :
    bo :
    supress_h :
    color : tuple or list of ints, optional, default None
        If tuple of ints all bonds selected will be highlighted with that color
        If list of OEColors, the list needs to be the same length as the incoming molecules
    names :

    Returns
    -------

    """
    itf = oechem.OEInterface()
    # PageByPage = True

    ropts = oedepict.OEReportOptions(rows, cols)
    ropts.SetHeaderHeight(25)
    ropts.SetFooterHeight(25)
    ropts.SetCellGap(2)
    ropts.SetPageMargins(10)
    report = oedepict.OEReport(ropts)

    cellwidth, cellheight = report.GetCellWidth(), report.GetCellHeight()
    opts = oedepict.OE2DMolDisplayOptions(cellwidth, cellheight, oedepict.OEScale_AutoScale)
    oedepict.OESetup2DMolDisplayOptions(opts, itf)

    if align:
        if isinstance(align, str):
            ref_mol = oechem.OEGraphMol()
            oechem.OESmilesToMol(ref_mol, align)
        elif isinstance(align, (oechem.OEMol, oechem.OEMolBase, oechem.OEGraphMol)):
            ref_mol = align
        oedepict.OEPrepareDepiction(ref_mol)

    for i, mol in enumerate(molecules):
        cell = report.NewCell()
        if isinstance(mol, str):
            m = oechem.OEMol()
            oechem.OESmilesToMol(m, mol)
            mol = m
            if names is not None:
                mol.SetTitle(str(names[i]))
        mol_copy = oechem.OEMol(mol)
        oedepict.OEPrepareDepiction(mol_copy, False, True)
        # if bo:
        #     b.SetData('WibergBondOrder', bo[i])
        #     opts.SetBondPropertyFunctor(LabelWibergBondOrder())

        atom_bond_set = oechem.OEAtomBondSet()
        a1 = mol_copy.GetAtom(oechem.OEHasMapIdx(bond_map_idx[0]))
        a2 = mol_copy.GetAtom(oechem.OEHasMapIdx(bond_map_idx[1]))
        b = mol_copy.GetBond(a1, a2)
        if wbos:
            b.SetData('WibergBondOrder', wbos[i])
        opts.SetBondPropertyFunctor(fragmenter.chemi.LabelWibergBondOrder())
        atom_bond_set.AddAtom(a1)
        atom_bond_set.AddAtom(a2)
        atom_bond_set.AddBond(b)
        hstyle = oedepict.OEHighlightStyle_BallAndStick
        if colors is not None:
            hcolor = oechem.OEColor(*colors[i])
        else:
            hcolor = oechem.OEColor(oechem.OELightBlue)

        overlaps = oegraphsim.OEGetFPOverlap(ref_mol, mol_copy, oegraphsim.OEGetFPType(oegraphsim.OEFPType_Tree))
        oedepict.OEPrepareMultiAlignedDepiction(mol_copy, ref_mol, overlaps)
        disp = oedepict.OE2DMolDisplay(mol_copy, opts)
        oedepict.OEAddHighlighting(disp, hcolor, hstyle, atom_bond_set)

        oedepict.OERenderMolecule(cell, disp)
        oedepict.OEDrawCurvedBorder(cell, oedepict.OELightGreyPen, 10.0)

    oedepict.OEWriteReport(fname, report)

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
    if maxY:
        sorted_indices = np.argsort(Xs)[::-1]
    else:
        sorted_indices = np.argsort(Xs)
    start = sorted_list[0]
    start.append(sorted_indices[0])
    pareto_front = [start]
    for i, pair in enumerate(sorted_list[1:]):
        if maxY:
            if pair[1] >= pareto_front[-1][1]:
                pair.append(sorted_indices[i+1])
                pareto_front.append(pair)

        else:
            if pair[1] <= pareto_front[-1][1]:
                pair.append(sorted_indices[i+1])
                pareto_front.append(pair)

    '''Plotting process'''
    fig, ax = plt.subplots()
    # ax.scatter(Xs,Ys,edgecolors='black', c='white')

    pf_X = [pair[0] for pair in pareto_front]
    pf_Y = [pair[1] for pair in pareto_front]
    plt.plot(pf_X, pf_Y, color='black', linewidth=0.8, zorder=1)

    ts = [0.001, 0.005, 0.01, 0.03, 0.05, 0.07, 0.1]
    if indices is not None:
        for i, j in enumerate(indices):
            j = int(j)
            if j == parent_idx:
                ax.scatter(Xs[j], Ys[j], marker='o', edgecolors='red', c='red', s=60, zorder=2)
                ax.annotate(ts[i], (Xs[j], Ys[j]), textcoords='offset points', xytext=(5, 4), color='red', fontsize=7)
                ax.annotate('   (parent)', (Xs[j], Ys[j]), textcoords='offset points', xytext=(5, 4), color='red',
                            fontsize=7)
            else:
                ax.scatter(Xs[j], Ys[j], marker='o', edgecolors='black',c='red', s=60, zorder=4)
                ax.annotate(ts[i], (Xs[j], Ys[j]), textcoords='offset points', xytext=(5, 4), color='black', fontsize=7)
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
    #plt.close()
    return pareto_front


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="calculate OE WBO")
    parser.add_argument('-n', '--name', type=str, help='Molecule name with number of its state')
    args = parser.parse_args()
    name = args.name

    # Load relevant data files
    with open('{}/{}_oe_wbo_by_bond.json'.format(name, name), 'r') as f:
        combinatorial_results = json.load(f)
    with open('{}/{}_wbo_dists.json'.format(name, name), 'r') as f:
        benchmark_results = json.load(f)

    for bond in benchmark_results:
        if bond == 'provenance':
            continue
        bond_des = fragmenter.utils.deserialize_bond(bond)

        # Find parent
        for smiles in combinatorial_results[bond]:
            if 'parent' in smiles:
                parent_smiles = combinatorial_results[bond][smiles]['map_to_parent']
                parent_wbos = combinatorial_results[bond][smiles]['individual_confs']

        parent_mol = oechem.OEMol()
        oechem.OESmilesToMol(parent_mol, parent_smiles)

        b = get_bond(parent_mol, bond_des)
        nbrs_1_4 = get_bond_nbr(parent_mol, b)
        frags_1_4 = frags_with_nbrs(combinatorial_results[bond], nbrs_1_4)
        scores = []
        costs = []
        frags = []
        map_to_parents = []
        elf10_wbos = []
        for i, smiles in enumerate(frags_1_4):
            if 'parent' in smiles:
                parent_idx = i
            score = mmd_x_xsqred(x=parent_wbos, y=combinatorial_results[bond][smiles]['individual_confs'])
            scores.append(score)
            smile = smiles.split('_')[0]
            cost = n_heavy_atoms(smile)
            costs.append(cost**3)
            frags.append(smile)
            map_to_parents.append(combinatorial_results[bond][smiles]['map_to_parent'])
            elf10_wbos.append(combinatorial_results[bond][smiles]['elf_estimate'])

        # Add results from benchmark
        benchmark_frags = {}
        benchmark_scores = []
        benchmark_costs = []
        for param in benchmark_results[bond]:
            params = param.split('_')
            if 'path' not in params:
                continue
            else:
                params.remove('length')
            threshold = params[0]
            hueristic = params[1]
            rotors = params[2]
            if rotors == 'True':
                continue
            if len(params) > 3:
                f = params[3]
                if f == 'False':
                    continue
            if threshold not in benchmark_frags:
                benchmark_frags[threshold] = benchmark_results[bond][param]['frag']
                benchmark_scores.append(mmd_x_xsqred(parent_wbos, benchmark_results[bond][param]['wbo_dist']))
                benchmark_costs.append(n_heavy_atoms(benchmark_results[bond][param]['frag']) ** 3)
            else:
                print('This should not happen')

        scores.extend(benchmark_scores)
        costs.extend(benchmark_costs)
        indices = np.linspace(len(scores)-7, len(scores)-1, 7)
        front = plot_pareto_frontier(Xs=scores, Ys=costs, parent_idx=parent_idx, maxY=False, indices=indices,
                                     filename='{}/{}_bond_{}_{}_front.pdf'.format(name, name, bond_des[0], bond_des[1]))

        # Visualize fragments on front
        # Get colors
        norm = plt.Normalize(0, max(scores))
        colors = plt.cm.viridis_r(norm(scores))

        front_frags = []
        colors_for_front_frags = []
        wbos = []
        front_scores = []
        for p in front:
            if p[-1] > len(map_to_parents)-1:
                # This is from the added scores
                continue
            front_frags.append(map_to_parents[p[-1]])
            colors_for_front_frags.append(rbg_to_int(colors[p[-1]], alpha=255))
            wbos.append(elf10_wbos[p[-1]])
            front_scores.append(p[0])
        to_pdf(front_frags, bond_map_idx=bond_des, wbos=wbos, colors=colors_for_front_frags, align=front_frags[0],
               names=front_scores, fname='{}/{}_bond_{}_{}_front_frags.pdf'.format(name, name, bond_des[0], bond_des[1]))
