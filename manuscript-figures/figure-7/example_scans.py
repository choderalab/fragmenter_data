import json
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
import pandas as pd
import math

from fragmenter import chemi
from openeye import oedepict, oechem, oegraphsim

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

    for i, s  in enumerate(smiles):
        cell = report.NewCell()
        mol = oechem.OEMol()
        oechem.OESmilesToMol(mol, s)
        oedepict.OEPrepareDepiction(mol, False, True)

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

        bond_label = oedepict.OEHighlightLabel("{:.2f}".format((wbos[i])), hcolor)
        oedepict.OEAddLabel(disp, bond_label, atom_bond_set)
        oedepict.OERenderMolecule(cell, disp)
        #oedepict.OEDrawCurvedBorder(cell, oedepict.OELightGreyPen, 10.0)

    return(oedepict.OEWriteReport(fname, report))

# Generate QC scan, WBO and highlight molecule
colors = chemi._KELLYS_COLORS
angles = np.arange(-180, 195, 15)

with open('../../phenyl_benchmark/data/qcarchive_torsiondrives.json', 'r') as f:
        fgroups_td = json.load(f)

fgroups = ['nitro', 'phenylurea']
indices = [[1, 3, 5], [0, 3, 8, 7]]
for fgroup, index in zip(fgroups, indices):

    smiles = [fgroups_td[fgroup]['indices'][i] for i in index]
    wbos = [fgroups_td[fgroup]['elf10_am1_wbo'][i] for i in index]
    energies = [fgroups_td[fgroup]['energy'][i] for i in index]
    lowdin_wbos = [fgroups_td[fgroup]['lowdin_wbos'][i] for i in index]
    visualize_phenyls(smiles, rows=5, cols=3, wbos=wbos, bond_idx=(2,3), colors=colors,
                       fname='figures/{}_example_td_drive_mols.pdf'.format(fgroup))

    plt.figure()
    linewidth = 1.0
    for i, e in enumerate(energies):
        plt.plot(angles, e, color=colors[i], linewidth=linewidth)
        plt.plot(angles, e, '.', color=colors[i])
        plt.xlabel('angles (degree)')
        plt.ylabel('relative energy (kJ/mol)')
    plt.savefig('figures/{}_example_qc_drives.pdf'.format(fgroup))
    plt.close()

    plt.figure()
    for i, w in enumerate(lowdin_wbos):
        plt.plot(angles[1:], w, color=colors[i], linewidth=linewidth)
        plt.plot(angles[1:], w, '*', color=colors[i])
        plt.xlabel('angles (degree)')
        plt.ylabel('Lowdin-Wiberg bond order')
    plt.savefig('figures/{}_example_lowdin_wiberg.pdf'.format(fgroup))
    plt.close()

