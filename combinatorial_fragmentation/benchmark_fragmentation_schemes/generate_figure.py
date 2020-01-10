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


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='look at WBO dist from an AM1 torsion scan')
    parser.add_argument('-n', '--name', type=str, help='Molecule name')
   # parser.add_argument('-b', '--bond', nargs= '+', type=int)
    args = parser.parse_args()
    name = args.name
    #bond = tuple(args.bond)

    #print(bond)

    with open('{}/{}_wbo_dists_fixed.json'.format(name, name), 'r') as f:
        results = json.load(f)
    with open('{}/{}_wbo_scans_fixed.json'.format(name, name), 'r') as f:
        scan_results = json.load(f)
    with open('{}/{}_pfizer_wbo_dists.json'.format(name, name), 'r') as f:
        pfizer_results = json.load(f)

    #
    # torsion_scans = {}
    # for bond in results:
    #     if bond == 'provenance':
    #         continue
    #     print(bond)
    #     bond_des = fragmenter.utils.deserialize_bond(bond)
    #     torsion_scans[bond] = {'parent': {'wbos': []}, 'pfizer': {'wbos': []}, 'wbo_scheme':{'wbos':[]}}
    #     for frag_type in torsion_scans[bond]:
    #         print(frag_type)
    #         mol = oechem.OEMol()
    #         if frag_type == 'pfizer':
    #             smiles = pfizer_results[bond]['frag']
    #             oe_wbo = pfizer_results[bond]['elf10_wbo']
    #         elif frag_type == 'wbo_scheme':
    #             key = '0.03_path_length_False_None'
    #             try:
    #                 smiles = results[bond][key]['frag']
    #                 oe_wbo = results[bond][key]['elf10_wbo']
    #             except KeyError:
    #                 key = '0.03_path_length_False'
    #                 smiles = results[bond][key]['frag']
    #                 oe_wbo = results[bond][key]['elf10_wbo']
    #         else:
    #             smiles = results[bond][frag_type]['frag']
    #             oe_wbo = results[bond][frag_type]['elf10_wbo']
    #         torsion_scans[bond][frag_type]['frag'] = smiles
    #         torsion_scans[bond][frag_type]['elf10_wbo'] = oe_wbo
    #
    #         mol = oechem.OEMol()
    #         oechem.OESmilesToMol(mol, smiles)
    #         dih = fragmenter.torsions.find_torsion_around_bond(molecule=mol, bond=bond_des)
    #         conformers = fragmenter.chemi.generate_grid_conformers(mol, dihedrals=[dih], intervals=[15], strict_types=False, strict_stereo=False)
    #         for conf in conformers.GetConfs():
    #             mol_copy = oechem.OEMol(conf)
    #             oechem.OEAddExplicitHydrogens(mol_copy)
    #             if oequacpac.OEAssignPartialCharges(mol_copy, oequacpac.OECharges_AM1BCCSym):
    #                 bo = get_bond(mol=mol_copy, bond_tuple=bond_des)
    #                 wbo = bo.GetData('WibergBondOrder')
    #                 torsion_scans[bond][frag_type]['wbos'].append(wbo)
    # # save wbos
    # with open('{}/{}_wbo_scans.json'.format(name, name), 'w') as f:
    #     json.dump(torsion_scans, f, indent=2, sort_keys=True)

    # Generate distribution

    for bond in results:
        if bond == 'provenance':
            continue

        des_bond = fragmenter.utils.deserialize_bond(bond)

        plt.figure()
        sbn.kdeplot(results[bond]['parent']['wbo_dist'], shade=True, label='parent molecule', color=sbn.color_palette('colorblind')[0])
        sbn.distplot(results[bond]['parent']['wbo_dist'], rug=True, hist=False, color=sbn.color_palette('colorblind')[0])

        score = mmd_x_xsqred(results[bond]['parent']['wbo_dist'], pfizer_results[bond]['wbo_dist'])
        sbn.kdeplot(pfizer_results[bond]['wbo_dist'], shade=True, color=sbn.color_palette('colorblind')[1],
                    label='Pfizer; score: {}'.format(round(score, 3)))
        sbn.distplot(pfizer_results[bond]['wbo_dist'], rug=True, hist=False, color=sbn.color_palette('colorblind')[1])

        score = mmd_x_xsqred(results[bond]['parent']['wbo_dist'], results[bond]['0.03']['wbo_dist'])
        sbn.kdeplot(results[bond]['0.03']['wbo_dist'], shade=True, color=sbn.color_palette('colorblind')[2],
                    label='WBO scheme; score: {}'.format(round(score, 3)))
        sbn.distplot(results[bond]['0.03']['wbo_dist'], rug=True, hist=False, color=sbn.color_palette('colorblind')[2])
        #sbn.distplot(results['0.03']['wbo_dist'], hist=False, color=sbn.color_palette()[2])

        plt.legend()
        plt.xticks(fontsize=14)
        plt.xlim(0.54, 1.45)
        plt.yticks([])
        plt.xlabel('Wiberg Bond Order', fontsize=14)
        plt.tight_layout()
        plt.savefig('{}/{}_bond_{}_{}_wbo_dist_fixed.pdf'.format(name, name, des_bond[0], des_bond[1]))

        # combine both scan and omega wbos
        plt.figure()
        x_parent = results[bond]['parent']['wbo_dist'] + scan_results[bond]['parent']['wbos']
        sbn.kdeplot(x_parent, shade=True, color=sbn.color_palette('colorblind')[0], label='parent molecule')
        sbn.distplot(x_parent, rug=True, hist=False, color=sbn.color_palette('colorblind')[0])

        x = pfizer_results[bond]['wbo_dist'] + scan_results[bond]['pfizer']['wbos']
        score = mmd_x_xsqred(x_parent, x)
        sbn.kdeplot(x, shade=True, color=sbn.color_palette('colorblind')[1], label='Pfizer; score: {}'.format(round(score, 3)))
        sbn.distplot(x, rug=True, hist=False, color=sbn.color_palette('colorblind')[1])

        x = results[bond]['0.03']['wbo_dist'] + scan_results[bond]['wbo_scheme']['wbos']
        score = mmd_x_xsqred(x_parent, x)
        sbn.kdeplot(x, shade=True, color=sbn.color_palette('colorblind')[2], label='WBO scheme; score: {}'.format(round(score, 3)))
        sbn.distplot(x, rug=True, hist=False, color=sbn.color_palette('colorblind')[2])

        plt.legend()
        plt.xticks(fontsize=14)
        plt.xlim(0.54)
        plt.yticks([])
        plt.xlabel('Wiberg Bond Order', fontsize=14)
        plt.tight_layout()
        print('figure')
        plt.savefig('{}/{}_bond_{}_{}_wbo_combined.pdf'.format(name, name, des_bond[0], des_bond[1]))

        smiles = [results[bond]['parent']['frag'], pfizer_results[bond]['frag'], results[bond]['0.03']['frag']]
        wbos = [results[bond]['parent']['elf10_wbo'], pfizer_results[bond]['elf10_wbo'], results[bond]['0.03']['elf10_wbo']]
        colors = [rbg_to_int(list(i), alpha=255) for i in sbn.color_palette('colorblind')[:3]]
        #colors.append(rbg_to_int(list(sbn.color_palette('colorblind')[4]), alpha=255))
        visualize_mols(smiles, cols=2, rows=2, bond_idx=des_bond, colors=colors, wbos=wbos,
                       fname='{}/{}_bond_{}_{}_frags_fixed_test.pdf'.format(name, name, des_bond[0], des_bond[1]),
                       align_to=2)


