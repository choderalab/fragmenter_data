import fragmenter
import json
from openeye import oechem, oequacpac, oedepict, oegraphsim
import matplotlib.pyplot as plt
import glob
import seaborn as sbn
import oenotebook as onb
import cmiles
import itertools
import numpy as np

def get_bond(mol, bond_tuple):
    a1 = mol.GetAtom(oechem.OEHasMapIdx(bond_tuple[0]))
    a2 = mol.GetAtom(oechem.OEHasMapIdx(bond_tuple[1]))
    if not a1 or not a2:
        print('no atoms')
        return False
    bond = mol.GetBond(a1, a2)
    if not bond:
        print('no bond')
        return False
    return bond

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
    minscale = float("inf")
    for s in smiles:
        mol = oechem.OEMol()
        oechem.OESmilesToMol(mol, s)
        mols.append(mol)
        oedepict.OEPrepareDepiction(mol, False, True)
        minscale = min(minscale, oedepict.OEGetMoleculeScale(mol, opts))
        print(minscale)

    print(minscale)
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
            hcolor = oechem.OEColor(*colors[i])

        overlaps = oegraphsim.OEGetFPOverlap(ref_mol, mol, oegraphsim.OEGetFPType(oegraphsim.OEFPType_Tree))
        oedepict.OEPrepareMultiAlignedDepiction(mol, ref_mol, overlaps)

        disp = oedepict.OE2DMolDisplay(mol, opts)
        oedepict.OEAddHighlighting(disp, hcolor, hstyle, atom_bond_set)

        bond_label = oedepict.OEHighlightLabel("{:.2f}".format((wbos[i])), hcolor)
        oedepict.OEAddLabel(disp, bond_label, atom_bond_set)
        oedepict.OERenderMolecule(cell, disp)
        # oedepict.OEDrawCurvedBorder(cell, oedepict.OELightGreyPen, 10.0)

    return (oedepict.OEWriteReport(fname, report))

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
name = 'Pemetrexed_0'
bond = '[7, 13]'
with open('{}_wbo_dists.json'.format(name), 'r') as f:
    results = json.load(f)
results = results[bond]
with open('{}_pfizer_wbo_dists.json'.format(name), 'r') as f:
    pfizer_results = json.load(f)

sbn.kdeplot(results['parent']['wbo_dist'], shade=True)
sbn.distplot(results['parent']['wbo_dist'], rug=True, hist=False, color=sbn.color_palette()[0])
sbn.distplot(results['parent']['wbo_dist'], hist=False, color=sbn.color_palette()[0])

sbn.kdeplot(results['0.01_path_length_False_None']['wbo_dist'], shade=True)
sbn.distplot(results['0.01_path_length_False_None']['wbo_dist'], rug=True, hist=False, color=sbn.color_palette()[1])
sbn.distplot(results['0.01_path_length_False_None']['wbo_dist'], hist=False, color=sbn.color_palette()[1])

sbn.kdeplot(results['0.03_path_length_False_None']['wbo_dist'], shade=True)
sbn.distplot(results['0.03_path_length_False_None']['wbo_dist'], rug=True, hist=False, color=sbn.color_palette()[2])
sbn.distplot(results['0.03_path_length_False_None']['wbo_dist'], hist=False, color=sbn.color_palette()[2])

sbn.kdeplot(results['0.1_path_length_False_None']['wbo_dist'], shade=True)
sbn.distplot(results['0.1_path_length_False_None']['wbo_dist'], rug=True, hist=False, color=sbn.color_palette()[3])
sbn.distplot(results['0.1_path_length_False_None']['wbo_dist'], hist=False, color=sbn.color_palette()[3])

sbn.kdeplot(pfizer_results[bond]['wbo_dist'], shade=True)
sbn.distplot(pfizer_results[bond]['wbo_dist'], rug=True, hist=False, color=sbn.color_palette()[4])
sbn.distplot(pfizer_results[bond]['wbo_dist'], hist=False, color=sbn.color_palette()[4])
plt.xticks(fontsize=14)
plt.yticks([])
plt.xlabel('Wiberg Bond Order', fontsize=14)
plt.tight_layout()
plt.savefig('wbo_dists.pdf')

colors = [rbg_to_int(list(i), alpha=255) for i in sbn.color_palette()]
wbos = [results['parent']['elf10_wbo'], results['0.01_path_length_False_None']['elf10_wbo'], results['0.03_path_length_False_None']['elf10_wbo'],
        results['0.1_path_length_False_None']['elf10_wbo'], pfizer_results[bond]['elf10_wbo']]
frags = [results['parent']['frag'], results['0.01_path_length_False_None']['frag'], results['0.03_path_length_False_None']['frag'],
         results['0.1_path_length_False_None']['frag'], pfizer_results[bond]['frag']]
des_bond = fragmenter.utils.deserialize_bond(bond)
visualize_mols(frags, cols=2, rows=2, bond_idx=des_bond, colors=colors, wbos=wbos, fname='fragments.pdf', align_to=0)