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
import qcfractal.interface as ptl
import cmiles

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

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)


def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


def find_improper_angles(mol):
    """
    Find the improper dihedral angles in some molecule. Currently supports
    those with a central trivalent nitrogen atom.
    Parameters
    ----------
    mol : OpenEye oemol
        oemol in which to look for improper angles
    Returns
    -------
    list
        Each element in the list is a 4-tuple of the coordinates for the
        atoms involved in the improper. The central atom is listed first
        in the tuple. Each member of the tuple is a numpy array.
    list
        List of strings for atoms in the improper, central atom is first.
    """
    mol_coords = mol.GetCoords()
    crdlist = []
    Idxlist = []
    for atom in mol.GetAtoms(oechem.OEIsInvertibleNitrogen()):
        # central atom
        aidx = atom.GetIdx()
        crd0 = np.asarray(mol_coords[aidx])
        # sort the neighbors
        nbors = sorted(list(atom.GetAtoms()))
        #check if there are 3 atoms connected to central atom in improper
        if len(nbors) != 3:
            return crdlist, namelist
        crd1 = np.asarray(mol_coords[nbors[0].GetIdx()])
        crd2 = np.asarray(mol_coords[nbors[1].GetIdx()])
        crd3 = np.asarray(mol_coords[nbors[2].GetIdx()])
        # store coordinates
        crdlist.append([crd0, crd1, crd2, crd3])
        Idxlist.append([atom.GetIdx(), nbors[0].GetIdx(), nbors[1].GetIdx(),nbors[2].GetIdx()])
    return crdlist, Idxlist


def calc_improper_angle(atom0, atom1, atom2, atom3, translate=True):
    """
    Calculate the improper dihedral angle of a set of given four atoms.
    Parameters
    ----------
    atom0 : numpy array
        CENTRAL atom coordinates
    atom1 : numpy array
        outer atom coordinates
    atom2 : numpy array
        outer atom coordinates
    atom3 : numpy array
        outer atom coordinates
    translate : bool
        True to translate central atom to origin, False to keep as is.
        This should not affect the results of the calculation.
    Returns
    -------
    float
        Angle in degrees.
    """
    if translate:
        atom1 = atom1 - atom0
        atom2 = atom2 - atom0
        atom3 = atom3 - atom0
        atom0 = atom0 - atom0 # central must be moved last
    # calculate vectors
    v0 = atom0-atom1
    v1 = atom2-atom1
    v2 = atom2-atom3
    w1 = np.cross(v0, v1)
    w2 = np.cross(v1, v2)
    angle = angle_between(w1,w2) # this angle should be in range [0,90]
    # compute distance from plane to central atom
    # eq 6 from http://mathworld.wolfram.com/Point-PlaneDistance.html
    # here I'm using atom1 for (x,y,z), but could also use atom2 or atom3
    numer = w2[0]*(atom0[0]-atom1[0]) + w2[1]*(atom0[1]-atom1[1]) + w2[2]*(atom0[2]-atom1[2])
    denom = np.sqrt(w2[0]**2 + w2[1]**2 + w2[2]**2)
    dist = numer/denom
    # set reference so that if central atom is above plane, angle -> [90,180]
    #print("this is the distance:" + str(dist))
    if dist < 0:
        angle = angle*-1
    #angle = abs(angle)
    return angle


# Generate QC scan, WBO and highlight molecule
colors = chemi._KELLYS_COLORS[3:]
angles = np.arange(-180, 195, 15)
symbols = ['o', 'P', '^', '*', 's', 'p',  'X', 'd', 'H', '>']

with open('../../phenyl_benchmark/data/qcarchive_torsiondrives.json', 'r') as f:
        fgroups_td = json.load(f)

fgroups = ['nitro', 'phenylurea', 'amino', 'carbamate']
indices = [[1, 3, 5], [0, 3, 8, 7], [0, 2, 8], [0, 2, 5]]
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
        plt.xlabel('Torsion angles (degree)')
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


# Plot pyramidal nitrogen for amino group
with open('../../phenyl_benchmark/data/qcarchive_torsiondrive_conformers.json', 'r') as f:
    td_conformers = json.load(f)

client = ptl.FractalClient()
ds = client.get_collection('TorsionDriveDataset', 'OpenFF Substituted Phenyl Set 1')

fgroup = 'amino'
index = [0, 2, 8]
smiles = [fgroups_td[fgroup]['indices'][i] for i in index]

all_angles = []
for sm in smiles:
    print(sm)
    a = []
    oemols = [cmiles.utils.mol_from_json(c) for c in td_conformers['amino'][sm]]
    entry = ds.get_entry(sm)
    dih = entry.td_keywords.dihedrals[0]
    print(dih)
    crds, indxs = find_improper_angles(oemols[0])

    if not crds:
        print('no trivalent nitrogen in {}'.format(fgroup))
        continue
    idx_to_use = None
    for j, idx in enumerate(indxs):
        n = 0
        for k in idx:
            if k in dih:
                n += 1
        if n == 3:
            idx_to_use = idx
            print(idx, dih)
    if idx_to_use == None:
        print('trivalent nitrogen is not in dihedral {}'.format(fgroup))
        continue
    for mol in td_conformers['amino'][sm]:
        coords = mol['geometry']
        coords = np.array(coords, dtype=float).reshape(int(len(coords)/3), 3)
        crds = [coords[i] for i in idx_to_use]
        #crds, idx = find_improper_angles(mol)
        #crds = crds[idx_to_use]
        angle = calc_improper_angle(crds[0], crds[1], crds[2], crds[3])
        a.append(angle)
    all_angles.append(a)

plt.figure()
angles = np.arange(-165, 195, 15)
for i, angle in enumerate(all_angles):
    plt.plot(angles, angle, symbols[i] ,color=colors[i])
    plt.plot(angles, angle ,color=colors[i])
    #plt.plot(angles, angles_0,  linewidth=0.8)
    plt.title('Improper angles along QC torsion scan')
    plt.xlabel('Torsion angle (degrees)')
    plt.ylabel('Improper angles (degrees)')
    #plt.xticks(fontsize=14)
    #plt.yticks(fontsize=14);
plt.tight_layout()
plt.savefig("example-improper-angles.pdf")