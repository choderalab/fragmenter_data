import fragmenter
import json
from openeye import oechem, oequacpac, oedepict, oegraphsim
import matplotlib.pyplot as plt
import glob
import seaborn as sbn
import cmiles
import itertools
import numpy as np

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

        #opts.SetBondPropLabelFontScale(4.0)
        disp = oedepict.OE2DMolDisplay(mol, opts)
        oedepict.OEAddHighlighting(disp, hcolor, hstyle, atom_bond_set)

        #font = oedepict.OEFont(oedepict.OEFontFamily_Default, oedepict.OEFontStyle_Bold, 12,
        #                       oedepict.OEAlignment_Default, oechem.OEBlack)
        bond_label = oedepict.OEHighlightLabel("{:.2f}".format((wbos[i])), hcolor)
        bond_label.SetFontScale(1.4)
        #bond_label.SetFont(font)

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

name = 'Ademetionine_0'
bonds = [(19, 10), (15, 6)]
optimal_smiles = ['[H][C@]1([C@H]([C@@H](O[C@@H]1C[S@@+](C)C)n2c3ncncc3nc2)O)[H]',
                  'C[S@@+](CC[C@H]([NH3+])C(=O)[O-])C']
optimal_mapped_smiles = ["[H:28][c:1]1[n:16][c:4]2[n:19]([c:2]([n:18][c:3]2[cH:5][n:17]1)[H:29])[C@@:10]3([O:24][C@@:9]([C@@:7]([C@:8]3([O:26][H:50])[H:31])([H:30])[H])([C:12]([S@+:27]([C:11]([H:34])([H:35])[H:36])[CH:14]([H:41])[H:42])([H:37])[H:38])[H:32])[H:33]",
    "[H:37][CH:12]([H:38])[S@+:27]([C:11]([H:34])([H:35])[H:36])[C:14]([H:41])([H:42])[C:13]([H:39])([H:40])[C@@:15]([H:43])([C:6](=[O:23])[O-:22])[N+:21]([H:46])([H:47])[H:48]"]


with open('{}_wbo_scans_fixed.json'.format(name), 'r') as f:
    torsion_scans = json.load(f)

with open('../../../combinatorial_fragmentation/rank_fragments/selected/{}/{}_oe_wbo_with_score.json'.format(name, name), 'r') as f:
    omega_results = json.load(f)
with open('{}_wbo_dists_fixed.json'.format(name), 'r') as f:
    omega_benchmark_results = json.load(f)
with open("{}_pfizer_wbo_dists.json".format(name), 'r') as f:
    pfizer_results = json.load(f)

for i, bond in enumerate(bonds):
    print(bond)
    ser_bond = fragmenter.utils.serialize_bond(bond)
    # Generate torsion scan for optimal fragment
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, optimal_mapped_smiles[i])
    dih = fragmenter.torsions.find_torsion_around_bond(mol, bond)
    torsion_scans[ser_bond]['optimal'] = {'wbos': [], 'elf10_wbo': omega_results[ser_bond][optimal_smiles[i]]['elf_estimate'],
                                      'frag': optimal_mapped_smiles[i]}
    conformers = fragmenter.chemi.generate_grid_conformers(mol, dihedrals=[dih], intervals=[15], strict_types=False)
    for conf in conformers.GetConfs():
        mol_copy = oechem.OEMol(conf)
        oechem.OEAddExplicitHydrogens(mol_copy)
        if oequacpac.OEAssignPartialCharges(mol_copy, oequacpac.OECharges_AM1BCCNoSym):
            bo = get_bond(mol=mol_copy, bond_tuple=bond)
            wbo = bo.GetData('WibergBondOrder')
            torsion_scans[ser_bond]['optimal']['wbos'].append(wbo)

    plt.figure()
    sbn.kdeplot(omega_benchmark_results[ser_bond]['parent']['wbo_dist'], shade=True, label='parent molecule')
    sbn.distplot(omega_benchmark_results[ser_bond]['parent']['wbo_dist'], rug=True, hist=False, color=sbn.color_palette('colorblind')[0])
    # sbn.distplot(results['parent']['wbo_dist'], hist=False, color=sbn.color_palette()[0])

    score = mmd_x_xsqred(omega_benchmark_results[ser_bond]['parent']['wbo_dist'], pfizer_results[ser_bond]['wbo_dist'])
    sbn.kdeplot(pfizer_results[ser_bond]['wbo_dist'], shade=True, label='Pfizer; score: {}'.format(round(score, 3)))
    sbn.distplot(pfizer_results[ser_bond]['wbo_dist'], rug=True, hist=False, color=sbn.color_palette()[1])
    # sbn.distplot(pfizer_results[bond]['wbo_dist'], hist=False, color=sbn.color_palette()[1])

    score = mmd_x_xsqred(omega_benchmark_results[ser_bond]['parent']['wbo_dist'], omega_benchmark_results[ser_bond]['0.03']['wbo_dist'])
    sbn.kdeplot(omega_benchmark_results[ser_bond]['0.03']['wbo_dist'], shade=True, label='WBO scheme; score: {}'.format(round(score, 3)))
    sbn.distplot(omega_benchmark_results[ser_bond]['0.03']['wbo_dist'], rug=True, hist=False, color=sbn.color_palette()[2])
    # sbn.distplot(results['0.03']['wbo_dist'], hist=False, color=sbn.color_palette()[2])

    score = mmd_x_xsqred(omega_benchmark_results[ser_bond]['parent']['wbo_dist'], omega_results[ser_bond][optimal_smiles[i]]['individual_confs'])
    sbn.kdeplot(omega_results[ser_bond][optimal_smiles[i]]['individual_confs'], shade=True,
                label='Optimal molecule; score: {}'.format(round(score, 3)))
    sbn.distplot(omega_results[ser_bond][optimal_smiles[i]]['individual_confs'], rug=True, hist=False, color=sbn.color_palette()[3])

    plt.legend()
    plt.xticks(fontsize=14)
    plt.xlim(0.54, 1.45)
    plt.yticks([])
    plt.xlabel('Wiberg Bond Order', fontsize=14)
    plt.tight_layout()
    plt.savefig('{}_bond_{}_{}_wbo_dist_optimal_fixed.pdf'.format( name, bond[0], bond[1]))

    # combine both scan and omega wbos
    plt.figure()
    x_parent = omega_benchmark_results[ser_bond]['parent']['wbo_dist'] + torsion_scans[ser_bond]['parent']['wbos']
    sbn.kdeplot(x_parent, shade=True, label='parent molecule')
    sbn.distplot(x_parent, rug=True, hist=False, color=sbn.color_palette()[0])

    x = pfizer_results[ser_bond]['wbo_dist'] + torsion_scans[ser_bond]['pfizer']['wbos']
    score = mmd_x_xsqred(x_parent, x)
    sbn.kdeplot(x, shade=True, label='Pfizer; score: {}'.format(round(score, 3)))
    sbn.distplot(x, rug=True, hist=False, color=sbn.color_palette()[1])

    x = omega_benchmark_results[ser_bond]['0.03']['wbo_dist'] + torsion_scans[ser_bond]['wbo_scheme']['wbos']
    score = mmd_x_xsqred(x_parent, x)
    sbn.kdeplot(x, shade=True, label='WBO scheme; score: {}'.format(round(score, 3)))
    sbn.distplot(x, rug=True, hist=False, color=sbn.color_palette()[2])

    x = omega_results[ser_bond][optimal_smiles[i]]['individual_confs'] + torsion_scans[ser_bond]['optimal']['wbos']

    score = mmd_x_xsqred(x_parent, x)
    sbn.kdeplot(x, shade=True,
                label='Optimal molecule; score: {}'.format(round(score, 3)))
    sbn.distplot(x, rug=True, hist=False,
                 color=sbn.color_palette()[3])

    plt.legend()
    plt.xticks(fontsize=14)
    plt.xlim(0.50, 1.45)
    plt.yticks([])
    plt.xlabel('Wiberg Bond Order', fontsize=14)
    plt.tight_layout()
    print('figure')
    plt.savefig('{}_bond_{}_{}_wbo_combined_optimal.pdf'.format( name, bond[0], bond[1]))

    smiles = [omega_benchmark_results[ser_bond]['parent']['frag'], pfizer_results[ser_bond]['frag'],
              omega_benchmark_results[ser_bond]['0.03']['frag'], optimal_mapped_smiles[i]]
    wbos = [omega_benchmark_results[ser_bond]['parent']['elf10_wbo'], pfizer_results[ser_bond]['elf10_wbo'],
            omega_benchmark_results[ser_bond]['0.03']['elf10_wbo'], omega_results[ser_bond][optimal_smiles[i]]['elf_estimate']]
    colors = [rbg_to_int(list(i), alpha=255) for i in sbn.color_palette()[:3]]
    visualize_mols(smiles, cols=2, rows=2, bond_idx=bond, colors=colors, wbos=wbos,
                   fname='{}_bond_{}_{}_frags_fixed_optimal.pdf'.format( name, bond[0], bond[1]),
                   align_to=2)
with open('{}_wbo_scans_optimal_fixed.json'.format(name), 'w') as f:
    json.dump(torsion_scans, f, indent=2, sort_keys=True)
