#!/usr/bin/env python
# coding: utf-8

import json
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
from matplotlib.ticker import FormatStrFormatter

import pandas as pd
import math
import arch.bootstrap
import seaborn as sbn

from fragmenter import chemi
from openeye import oedepict, oechem, oegraphsim

sbn.set_style('whitegrid')
sbn.set_context('paper', font_scale=1.7)

def get_bond(mol, bond_idx):
    atoms = [mol.GetAtom(oechem.OEHasMapIdx(i)) for i in bond_idx]
    bond = mol.GetBond(atoms[0], atoms[1])
    return bond

def visualize_phenyls(smiles, fname, rows, cols, bond_idx, wbos, colors):
    """
    Visualize molecules with highlighted bond and labeled with WBO
    Parameters
    ----------
    smiles : list of SMILES to visualize. Torsion atoms should have map indices (1, 2, 3, 4)
    fname : str
        filename
    rows : int
    cols : int
    bond_idx : tuple of atom maps of bond to highlight. Since all torsions are mapped with (1, 2, 3, 4), the same tuple
    is used for all SMILES.
    wbos : list of floats
    colors : list of hex values for colors

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

    pen = oedepict.OEPen(oechem.OEBlack, oechem.OEBlack, oedepict.OEFill_Off, 0.9)
    opts.SetDefaultBondPen(pen)

    # align to first molecule
    ref_mol = oechem.OEGraphMol()
    oechem.OESmilesToMol(ref_mol, smiles[-1])
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
        #mol = oechem.OEMol()
        #oechem.OESmilesToMol(mol, s)
        oedepict.OEPrepareDepiction(mol, False, True)

        bond = get_bond(mol, bond_idx)
        atom_bond_set = oechem.OEAtomBondSet()
        atom_bond_set.AddAtoms([bond.GetBgn(), bond.GetEnd()])
        atom_bond_set.AddBond(bond)

        hstyle = oedepict.OEHighlightStyle_BallAndStick
        hcolor = oechem.OEColor(colors[len(mols)-i-1])

        overlaps = oegraphsim.OEGetFPOverlap(ref_mol, mol, oegraphsim.OEGetFPType(oegraphsim.OEFPType_Tree))
        oedepict.OEPrepareMultiAlignedDepiction(mol, ref_mol, overlaps)

        disp = oedepict.OE2DMolDisplay(mol, opts)
        oedepict.OEAddHighlighting(disp, hcolor, hstyle, atom_bond_set)

        font = oedepict.OEFont(oedepict.OEFontFamily_Default, oedepict.OEFontStyle_Default, 24,
                               oedepict.OEAlignment_Default, oechem.OEBlack)
        bond_label = oedepict.OEHighlightLabel("{:.2f}".format(wbos[i]), hcolor)
        bond_label.SetFontScale(4.0)
        bond_label.SetFont(font)
        oedepict.OEAddLabel(disp, bond_label, atom_bond_set)
        oedepict.OERenderMolecule(cell, disp)
        #oedepict.OEDrawCurvedBorder(cell, oedepict.OELightGreyPen, 10.0)

    return(oedepict.OEWriteReport(fname, report))


with open('data/qcarchive_torsiondrives.json', 'r') as f:
    fgroups_td = json.load(f)


fgroups_smarts = {
    #'phenoxide': 'C[O-]',
    'dimethylamino': 'N(C)(C)',
    'methylamino': 'NC',
    'amino': 'N',
    'ethylamino': 'NCC',
    'propylamino': 'NCCC',
    'hydroxy': 'O',
    'methoxy': 'OC',
    'ethoxy': 'OCC',
    'dimethylurea': 'NC(=O)N(C)(C)',
    'urea': 'NC(=O)NC',
    'phenylurea': 'NC(=O)N',
    'ethylamide': 'NC(=O)CC',
    'amide': 'NC(=O)C',
    #'fluoro': 'CF',
    #'chloro': 'CCl',
    #'methyl': 'C',
    #'cyano': 'CC#N',
    #'bromo': 'CBr',
    'carbamate': 'OC(=O)N',
    #'iodo': 'CI',
    'benzoicacid': 'C(=O)O',
    'ethoxycarbonyl': 'C(=O)OCC',
    #'trifluoromethyl': 'CC(F)(F)(F)',
    #'trimethylamonium': 'C[N+](C)(C)C',
    'nitro': '[N+](=O)[O-]'
}
color_keys = ['rosybrown', 'indianred', 'red', 'orange', 'gold', 'yellow','greenyellow', 'green', 'limegreen',
          'lightseagreen', 'teal', 'cyan', 'deepskyblue', 'mediumslateblue', 'blueviolet', 'mediumorchid', 'lightpink']
colors = mcolors.CSS4_COLORS


# Generate regression plot of WBO and torsion energy barrier height
for i, fgroup in enumerate(fgroups_smarts):
    energies = fgroups_td[fgroup]['energy']
    am1_wbos = fgroups_td[fgroup]['elf10_am1_wbo']
    max_energies = [max(energy) for energy in energies]
    slope, intercept, r_value, p_value, std_err = stats.linregress(am1_wbos, max_energies)
    fgroups_td[fgroup]['stats'] = [slope, r_value**2, p_value, std_err]
    plt.plot(np.unique(am1_wbos), np.poly1d([slope, intercept])(np.unique(am1_wbos)), color_keys[i])
    plt.scatter(x=am1_wbos, y=max_energies, color=color_keys[i], label=fgroups_smarts[fgroup], s=4)
plt.legend(bbox_to_anchor=(1, 1))
plt.xlabel('ELF10 AM1 Wiberg bond order')
plt.ylabel('Torsion barrier height (kJ/mol)')
plt.tight_layout()
plt.savefig('figures/qcarchive_torsiondrives/energy_vs_wbo.pdf')


# generate table with regression statistics
stats_table = {'functional group': [], 'slope': [], 'r^2': [], 'P value': [], 'standard error': []}
for fgroup in fgroups_smarts:
    stats_table['functional group'].append(fgroups_smarts[fgroup])
    stats_table['slope'].append(fgroups_td[fgroup]['stats'][0])
    stats_table['r^2'].append(fgroups_td[fgroup]['stats'][1])
    stats_table['P value'].append(fgroups_td[fgroup]['stats'][2])
    stats_table['standard error'].append(fgroups_td[fgroup]['stats'][3])
latex_table = pd.DataFrame(stats_table).to_latex(index=False)
with open('figures/qcarchive_torsiondrives/stats.tex', 'w') as f:
    f.write(latex_table)

# Generate QC scan, WBO and highlight molecule
colors = chemi._KELLYS_COLORS
angles = np.arange(-180, 195, 15)

def r_value_ci(am1_wbos, max_energies):
    return stats.linregress(am1_wbos, max_energies)[2]**2

for j, fgroup in enumerate(fgroups_smarts):
    print(fgroup)
    n = len(fgroups_td[fgroup]['indices'])
    cols = 3
    rows = math.ceil(n/cols)
    visualize_phenyls(list(reversed(fgroups_td[fgroup]['indices'])), rows=6, cols=2,
                        wbos=list(reversed(fgroups_td[fgroup]['elf10_am1_wbo'])), bond_idx=(2,3), colors=colors,
                       fname='figures/qcarchive_torsiondrives/{}_driven_mols.pdf'.format(fgroup))
     
    fig, ax = plt.subplots()
    #ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    linewidth = 1.0
    for i, e in enumerate(fgroups_td[fgroup]['energy']):
        plt.plot(angles, e, color=colors[i], linewidth=linewidth)
        plt.plot(angles, e, '.', color=colors[i])
        plt.xlabel('angles (degree)')
        plt.ylabel('relative energy (kJ/mol)')
        plt.xlim(-180, 180)
    plt.title('QC torsion scans')
    fig.tight_layout()
    plt.savefig('figures/qcarchive_torsiondrives/{}_qc_drives.pdf'.format(fgroup), bbox='tight')
    plt.close()

    fig, ax = plt.subplots()
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    for i, w in enumerate(fgroups_td[fgroup]['lowdin_wbos']):
        plt.plot(angles[1:], w, color=colors[i], linewidth=linewidth)
        plt.plot(angles[1:], w, '*', color=colors[i])
        plt.xlabel('angles (degree)')
        plt.ylabel('Wiberg bond order')
        plt.xlim(-180, 180)
    plt.title('WBO along torsion scans')
    fig.tight_layout()
    plt.savefig('figures/qcarchive_torsiondrives/{}_lowdin_wiberg.pdf'.format(fgroup), bbbox='tight')
    plt.close()
    
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    energies = fgroups_td[fgroup]['energy']
    max_energies = [max(energy) for energy in energies]
    am1_wbos = fgroups_td[fgroup]['elf10_am1_wbo']
    slope, intercept, r_value, p_value, std_err = stats.linregress(am1_wbos, max_energies)
    r_ci = arch.bootstrap.IIDBootstrap(np.asarray(am1_wbos), np.asarray(max_energies)).conf_int(r_value_ci, 10000,
                                                                                                method='percentile')
    plt.plot(np.unique(am1_wbos), np.poly1d([slope, intercept])(np.unique(am1_wbos)), 'black')
    for i, (wbo, en) in enumerate(zip(am1_wbos, max_energies)):
        ax.scatter(wbo, en, c=colors[i], s=50)
    textstr = '\n'.join((
        r'slope=%.2f kJ/mol' % (slope, ),
        r'$r^2=%.2f_{%.2f}^{%.2f}$' % (r_value**2, r_ci[0][0], r_ci[1][0])))
        #r'P value=%.3f' % (p_value, ),
        #r'standard error=%.2f' % (std_err, )))
    props = dict(boxstyle='square', facecolor='white', alpha=0.5)
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, #fontsize=12,
            verticalalignment='top', bbox=props)
    plt.xlabel('ELF10 AM1 Wiberg Bond Order')
    plt.ylabel('Energy Barrier height (kJ/mol)')
    ax.set_title('WBO vs torsion energy barrier')
    fig.tight_layout()
    plt.savefig('figures/qcarchive_torsiondrives/{}_energy_vs_wbo.pdf'.format(fgroup), bbox='tight')
    plt.close()