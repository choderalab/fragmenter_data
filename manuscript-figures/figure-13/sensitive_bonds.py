import fragmenter
import json
from openeye import oechem,  oedepict
import matplotlib.pyplot as plt
import matplotlib as mpl
import glob
import numpy as np
import itertools

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

    mmd2 = np.sqrt((E_x - E_y)**2 + (E_x_squared - E_y_squared)**2)
    return mmd2

def get_bond(mol, bond_idx):
    a1 = mol.GetAtom(oechem.OEHasMapIdx(bond_idx[0]))
    a2 = mol.GetAtom(oechem.OEHasMapIdx(bond_idx[1]))
    bond = mol.GetBond(a1, a2)
    if not bond:
        print('bond {} not connnected'.format(bond_idx))
        return False
        #raise ValueError("({}) atoms are not connected".format(bond_idx))
    return bond


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


def visualize_bond_atom_sensitivity(mols, bonds, scores, fname, rows, cols, atoms=None, min_scale=True):
    """

    Parameters
    ----------
    mols :
    bonds :
    scores :
    fname :
    wbos :
    rows :
    cols :
    atoms :
    height :
    width :

    Returns
    -------

    """

    itf = oechem.OEInterface()
    ropts = oedepict.OEReportOptions(rows, cols)
    ropts.SetHeaderHeight(0.01)
    ropts.SetFooterHeight(0.01)
    ropts.SetCellGap(0.0001)
    ropts.SetPageMargins(0.01)
    report = oedepict.OEReport(ropts)

    cellwidth, cellheight = report.GetCellWidth(), report.GetCellHeight()
    opts = oedepict.OE2DMolDisplayOptions(cellwidth, cellheight, oedepict.OEScale_AutoScale)
    oedepict.OESetup2DMolDisplayOptions(opts, itf)
    opts.SetAromaticStyle(oedepict.OEAromaticStyle_Circle)
    opts.SetAtomColorStyle(oedepict.OEAtomColorStyle_WhiteMonochrome)

    pen = oedepict.OEPen(oechem.OEBlack, oechem.OEBlack, oedepict.OEFill_Off, 0.9)
    opts.SetDefaultBondPen(pen)
    oedepict.OESetup2DMolDisplayOptions(opts, itf)

    if min_scale:
        minscale = float("inf")
        for m in mols:
            oedepict.OEPrepareDepiction(m, False, True)
            minscale = min(minscale, oedepict.OEGetMoleculeScale(m, opts))

        opts.SetScale(minscale)
    for i, mol in enumerate(mols):
        cell = report.NewCell()
        oedepict.OEPrepareDepiction(mol, False, True)
        atom_bond_sets = []
        for j, bond in enumerate(bonds[i]):
            bo = get_bond(mol, bond)
            atom_bond_set = oechem.OEAtomBondSet()
            atom_bond_set.AddBond(bo)
            atom_bond_sets.append(atom_bond_set)

        opts.SetTitleLocation(oedepict.OETitleLocation_Hidden)

        disp = oedepict.OE2DMolDisplay(mol, opts)
        hstyle = oedepict.OEHighlightStyle_Stick
        hstyle_2 = oedepict.OEHighlightStyle_Color
        score = scores[i]
        norm = plt.Normalize(0, max(score))
        colors = plt.cm.coolwarm(norm(score))
        colors_oe = [rbg_to_int(c, 200) for c in colors]

        for j, atom_bond_set in enumerate(atom_bond_sets):
            highlight = oechem.OEColor(*colors_oe[j])

            oedepict.OEAddHighlighting(disp, highlight, hstyle, atom_bond_set)
            oedepict.OEAddHighlighting(disp, highlight, hstyle_2, atom_bond_set)

        highlight = oedepict.OEHighlightByCogwheel(oechem.OEDarkPurple)
        highlight.SetBallRadiusScale(5.0)

        if not atoms is None:
            for a_b in atoms[i]:
                if isinstance(a_b[-1], list):
                    for k, c in enumerate(a_b[-1]):
                        print(c)
                        color = oechem.OEColor(*colors_oe[c])
                        highlight.SetBallRadiusScale(5.0 - 2.5 * k)
                        highlight.SetColor(color)
                        atom_bond_set_a = oechem.OEAtomBondSet()
                        if len(a_b[0]) == 1:
                            a = mol.GetAtom(oechem.OEHasMapIdx(a_b[0][0]))
                            atom_bond_set_a.AddAtom(a)
                        oedepict.OEAddHighlighting(disp, highlight, atom_bond_set_a)
                else:
                    color = oechem.OEColor(*colors_oe[a_b[-1]])
                    highlight.SetColor(color)
                    atom_bond_set_a = oechem.OEAtomBondSet()
                    if len(a_b[0]) == 1:
                        a = mol.GetAtom(oechem.OEHasMapIdx(a_b[0][0]))
                        atom_bond_set_a.AddAtom(a)
                    else:
                        for b in itertools.combinations(a_b[0], 2):
                            bo = get_bond(mol, b)
                            if not bo:
                                continue
                            atom_bond_set_a.AddAtom(bo.GetBgn())
                            atom_bond_set_a.AddAtom(bo.GetEnd())
                            atom_bond_set_a.AddBond(bo)
                    oedepict.OEAddHighlighting(disp, highlight, atom_bond_set_a)
        oedepict.OERenderMolecule(cell, disp)
        # oedepict.OEDrawCurvedBorder(cell, oedepict.OELightGreyPen, 10.0)

    return oedepict.OEWriteReport(fname, report)

def generate_figure(mol_names, fname, rows=4, cols=3, min_scale=True):
    all_bonds = []
    all_scores = []
    lengths = []
    all_smiles = []
    names = []
    for name in mol_names:
        with open(
                '../../combinatorial_fragmentation/rank_fragments/selected/{}/{}_oe_wbo_with_score.json'.format(name, name),
                'r') as f:
            results = json.load(f)

        bonds = []
        scores = []
        for bond in results:
            bond_des = fragmenter.utils.deserialize_bond(bond)
            for key in results[bond]:
                if 'parent' in key:
                    parent_key = key
                    parent_wbos = results[bond][parent_key]['individual_confs']
                    parent_smiles = results[bond][parent_key]['map_to_parent']
                    continue
            score = 0
            for key in results[bond]:
                score = max(score, mmd_x_xsqred(parent_wbos, results[bond][key]['individual_confs']))
            bonds.append(bond_des)
            scores.append(score)
        all_bonds.append(bonds)
        all_scores.append(scores)
        all_smiles.append(parent_smiles)
        if len(lengths) == 0:
            lengths.append(len(bonds))
        else:
            lengths.append(len(bonds) + lengths[-1])
        names.append(name)
        mols = []
        for s in all_smiles:
            mol = oechem.OEMol()
            oechem.OESmilesToMol(mol, s)
            mols.append(mol)

    visualize_bond_atom_sensitivity(mols, bonds=all_bonds, scores=all_scores, rows=rows, cols=cols, min_scale=min_scale,
                                        fname=fname)
    return all_scores

# Generate SI with all molecules
mol_names = glob.glob('../../combinatorial_fragmentation/rank_fragments/selected/*')
all_mols = [name.split('/')[-1] for name in mol_names]
generate_figure(all_mols, rows=5, cols=5, min_scale=False, fname='all_selected_mols.pdf')

#ToDo Generate colorbars for molecules in panel 1 and 2.
panel_1_mols = [ 'Ademetionine_0', 'Almitrine_1', 'Amlodipine_0', 'Bosutinib_0', 'Ceftazidime_0', 'Eltrombopag_1',
                 'Sulfinpyrazone_0', 'Fostamatinib_0', 'Nizatidine_0']

all_scores = generate_figure(panel_1_mols, min_scale=True, rows=4, cols=3, fname='panel_1.pdf')
# Generate colorbars to use in figure
font_size = 14
fig, axs = plt.subplots(6, 3)
n = 0
for i in range(5):
    for j in range(3):
        print(n)
        if n > 8:
            continue
        #f = plt.subplot(1, i+1)
        cmap = mpl.cm.coolwarm
        norm = mpl.colors.Normalize(0, max(all_scores[n]))
        cb1 = mpl.colorbar.ColorbarBase(axs[i, j], cmap=cmap,
                            norm=norm,
                        orientation='horizontal')
        cb1.ax.tick_params(labelsize=font_size)
        cb1.set_label(n)
        mi = 0
        mx = max(all_scores[n])
        mid = (mi + mx)/2
        cb1.set_ticks([mi, mid, mx])
        cb1.set_ticklabels([mi, round(mid, 2), round(mx, 2)])
        #cb1.ax.locator_params(nbins=2)
        n+=1
plt.tight_layout()
plt.savefig('colorbars_panel_1.pdf')


panel_2_mols = ['Dacomitinib_0', 'Tedizolid_phosphate_0', 'Acemetacin_0']


all_bonds = []
all_scores = []
lengths = []
all_smiles = []
names = []
for name in panel_2_mols:
    with open(
            '../../combinatorial_fragmentation/rank_fragments/selected/{}/{}_oe_wbo_with_score.json'.format(name, name),
            'r') as f:
        results = json.load(f)

    bonds = []
    scores = []
    for bond in results:
        bond_des = fragmenter.utils.deserialize_bond(bond)
        for key in results[bond]:
            if 'parent' in key:
                parent_key = key
                parent_wbos = results[bond][parent_key]['individual_confs']
                parent_smiles = results[bond][parent_key]['map_to_parent']
                continue
        score = 0
        for key in results[bond]:
            score = max(score, mmd_x_xsqred(parent_wbos, results[bond][key]['individual_confs']))
        bonds.append(bond_des)
        scores.append(score)
    all_bonds.append(bonds)
    all_scores.append(scores)
    all_smiles.append(parent_smiles)
    if len(lengths) == 0:
        lengths.append(len(bonds))
    else:
        lengths.append(len(bonds) + lengths[-1])
    names.append(name)
    mols = []
    for s in all_smiles:
        mol = oechem.OEMol()
        oechem.OESmilesToMol(mol, s)
        mols.append(mol)

atoms = [[([29], all_bonds[0].index((28, 17))), ([30], all_bonds[0].index((27, 14)))],
             [([24, 25, 27, 29, 31], all_bonds[1].index((9, 23)))],
             [([23], [all_bonds[2].index((17, 28)), all_bonds[2].index((15, 22))])]]

visualize_bond_atom_sensitivity(mols, bonds=all_bonds, scores=all_scores, rows=4, cols=3, atoms=atoms, min_scale=True,
                                    fname='panel_2.pdf')

fig, axs = plt.subplots(6, 3)
n = 0
for i in range(5):
    for j in range(3):
        if n > 2:
            continue
        print(n)
        #f = plt.subplot(1, i+1)
        # f = plt.subplot(1, i+1)
        cmap = mpl.cm.coolwarm
        norm = mpl.colors.Normalize(0, max(all_scores[n]))
        cb1 = mpl.colorbar.ColorbarBase(axs[i, j], cmap=cmap,
                                        norm=norm,
                                        orientation='horizontal')
        cb1.ax.tick_params(labelsize=font_size)
        cb1.set_label(n)
        mi = 0
        mx = max(all_scores[n])
        mid = (mi + mx) / 2
        cb1.set_ticks([mi, mid, mx])
        cb1.set_ticklabels([mi, round(mid, 2), round(mx, 2)])
        # cb1.ax.locator_params(nbins=2)
        n += 1
plt.tight_layout()
plt.savefig('colorbars_panel_2.pdf')