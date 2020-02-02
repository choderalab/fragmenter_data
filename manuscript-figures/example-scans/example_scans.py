import json
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

import seaborn as sbn
import numpy as np
import matplotlib.colors as mcolors
import pandas as pd
import math

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

    # align to first molecule
    ref_mol = oechem.OEGraphMol()
    oechem.OESmilesToMol(ref_mol, smiles[0])
    oedepict.OEPrepareDepiction(ref_mol)

    mols = []
    minscale = float("inf")
    for s in smiles:
        mol = oechem.OEMol()
        oechem.OESmilesToMol(mol, s)
        mols.append(mol)
        oedepict.OEPrepareDepiction(mol, False, True)
        minscale = min(minscale, oedepict.OEGetMoleculeScale(mol, opts))
        opts.SetScale(minscale)

    print(minscale)
    opts.SetScale(minscale)

    for i, mol  in enumerate(mols):
        cell = report.NewCell()
        #mol = oechem.OEMol()
        #oechem.OESmilesToMol(mol, s)
        #oedepict.OEPrepareDepiction(mol, False, True)

        bond = get_bond(mol, bond_idx)
        atom_bond_set = oechem.OEAtomBondSet()
        atom_bond_set.AddAtoms([bond.GetBgn(), bond.GetEnd()])
        atom_bond_set.AddBond(bond)

        hstyle = oedepict.OEHighlightStyle_BallAndStick
        hcolor = oechem.OEColor(colors[i])

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

# Generate QC scan, WBO and highlight molecule
colors = chemi._KELLYS_COLORS[3:]
angles = np.arange(-180, 195, 15)
symbols = ['o', 'P', '^', '*', 's', 'p',  'X', 'd', 'H', '>']

with open('../../phenyl_benchmark/data/qcarchive_torsiondrives.json', 'r') as f:
        fgroups_td = json.load(f)

fgroups = ['nitro', 'phenylurea']
indices = [[1, 3, 5], [0, 3, 8, 7]]
for fgroup, index in zip(fgroups, indices):
    print(fgroup)

    smiles = [fgroups_td[fgroup]['indices'][i] for i in index]
    wbos = [fgroups_td[fgroup]['elf10_am1_wbo'][i] for i in index]
    energies = [fgroups_td[fgroup]['energy'][i] for i in index]
    lowdin_wbos = [fgroups_td[fgroup]['lowdin_wbos'][i] for i in index]
    # Add last to first (it's what I did for the energies so it goes from -180 - 180
    visualize_phenyls(smiles, rows=2, cols=2, wbos=wbos, bond_idx=(2,3), colors=colors,
                       fname='{}_example_td_drive_mols.pdf'.format(fgroup))

    plt.figure()
    linewidth = 1.0
    for i, e in enumerate(energies):
        plt.plot(angles, e, color=colors[i], linewidth=linewidth)
        plt.plot(angles, e, symbols[i], color=colors[i])
        plt.xlabel('angles (degree)')
        plt.ylabel('relative energy (kJ/mol)')
    plt.tight_layout()
    plt.ylim(-0.5, 40)
    plt.savefig('{}_example_qc_drives.pdf'.format(fgroup))
    plt.close()

    fig, ax = plt.subplots()
    for i, w in enumerate(lowdin_wbos):
        # Insert last to first so it coincides with what I did to the energy scan (this is so that the scan goes from -180-180)
        w.insert(0, w[-1])
        plt.plot(angles, w, color=colors[i], linewidth=linewidth)
        plt.plot(angles, w, symbols[i], color=colors[i])
        plt.xlabel('angles (degree)')
        plt.ylabel('Wiberg bond order')
    plt.tight_layout()
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.savefig('{}_example_lowdin_wiberg.pdf'.format(fgroup))
    plt.close()

    # Generate correlation plot
    fig, ax = plt.subplots()
    for i, values in enumerate(zip(energies, lowdin_wbos)):
        corr = np.corrcoef(values[1], values[0])[0][1]
        ax.scatter(values[1], values[0], marker=symbols[i], color=colors[i], label=r'$\rho =%.2f$' % (corr))
        ax.set_xlabel('Wiberg Bond Order')
        ax.set_ylabel('Conformer energy (kJ/mol)')

    l = ax.legend(loc='upper right', bbox_to_anchor=(1, 1), framealpha=0.5, fancybox=True)
    print(l)
    for i, text in enumerate(l.get_texts()):
        text.set_color(colors[i])
    plt.tight_layout()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    plt.savefig('{}_scatter.pdf'.format(fgroup))
    plt.close()

