import fragmenter
import json
from openeye import oechem, oequacpac, oedepict, oegraphsim
import matplotlib.pyplot as plt
import glob
import seaborn as sbn
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


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='look at WBO dist from an AM1 torsion scan')
    parser.add_argument('-n', '--name', type=str, help='Molecule name')
    args = parser.parse_args()
    name = args.name

    with open('{}/{}_wbo_dists_fixed.json'.format(name, name), 'r') as f:
        results = json.load(f)
    with open('{}/{}_pfizer_wbo_dists.json'.format(name, name), 'r') as f:
        pfizer_results = json.load(f)
    with open('{}/{}_wbo_scans.json'.format(name, name), 'r') as f:
        scan_results = json.load(f)

    torsion_scans = {}
    for bond in results:
        if bond == 'provenance':
            continue
        print(bond)
        bond_des = fragmenter.utils.deserialize_bond(bond)
        torsion_scans[bond] = {'parent': {'wbos': []}, 'pfizer': {'wbos': []}, 'wbo_scheme':{'wbos':[]}}
        for frag_type in torsion_scans[bond]:
            print(frag_type)
            mol = oechem.OEMol()
            if frag_type == 'pfizer':
                smiles = pfizer_results[bond]['frag']
                oe_wbo = pfizer_results[bond]['elf10_wbo']
            elif frag_type == 'wbo_scheme':
                key = '0.03'
                #try:
                smiles = results[bond][key]['frag']
                oe_wbo = results[bond][key]['elf10_wbo']
                # except KeyError:
                #     key = '0.03_path_length_False'
                #     smiles = results[bond][key]['frag']
                #     oe_wbo = results[bond][key]['elf10_wbo']
            else:
                smiles = results[bond][frag_type]['frag']
                oe_wbo = results[bond][frag_type]['elf10_wbo']
            torsion_scans[bond][frag_type]['frag'] = smiles
            torsion_scans[bond][frag_type]['elf10_wbo'] = oe_wbo

            if smiles == scan_results[bond][frag_type]['frag']:
                print('{} already scanned'.format(smiles))
                torsion_scans[bond][frag_type]['wbos'] = scan_results[bond][frag_type]['wbos']
                continue
            print('{} not found. scanning...'.format(smiles))
            mol = oechem.OEMol()
            oechem.OESmilesToMol(mol, smiles)
            dih = fragmenter.torsions.find_torsion_around_bond(molecule=mol, bond=bond_des)
            conformers = fragmenter.chemi.generate_grid_conformers(mol, dihedrals=[dih], intervals=[15], strict_types=False, strict_stereo=False)
            for conf in conformers.GetConfs():
                mol_copy = oechem.OEMol(conf)
                oechem.OEAddExplicitHydrogens(mol_copy)
                if oequacpac.OEAssignPartialCharges(mol_copy, oequacpac.OECharges_AM1BCCSym):
                    bo = get_bond(mol=mol_copy, bond_tuple=bond_des)
                    wbo = bo.GetData('WibergBondOrder')
                    torsion_scans[bond][frag_type]['wbos'].append(wbo)
    # save wbos
    with open('{}/{}_wbo_scans_fixed.json'.format(name, name), 'w') as f:
        json.dump(torsion_scans, f, indent=2, sort_keys=True)
    # with open('{}/{}_wbo_scans.json'.format(name, name), 'r') as f:
    #     torsion_scans = json.load(f)
    # Generate distribution
    for bond in torsion_scans:
        plt.figure()
        bond_des = fragmenter.utils.deserialize_bond(bond)
        for i, frag_type in enumerate(torsion_scans[bond]):
            sbn.kdeplot(torsion_scans[bond][frag_type]['wbos'], shade=True)
            sbn.distplot(torsion_scans[bond][frag_type]['wbos'], rug=True, hist=False, color=sbn.color_palette()[i])
            sbn.distplot(torsion_scans[bond][frag_type]['wbos'], hist=False, color=sbn.color_palette()[i])

        plt.xticks(fontsize=14)
        plt.xlim(0.54, 1.5)
        plt.yticks([])
        plt.xlabel('Wiberg Bond Order', fontsize=14)
        plt.tight_layout()
        plt.savefig('{}/{}_bond_{}_{}_wbo_scan_dist_fixed.pdf'.format(name, name, bond_des[0], bond_des[1]))

        colors = [rbg_to_int(list(i), alpha=255) for i in sbn.color_palette()[:3]]
        wbos = []
        frags = []
        for frag_type in torsion_scans[bond]:
            wbos.append(torsion_scans[bond][frag_type]['elf10_wbo'])
            frags.append(torsion_scans[bond][frag_type]['frag'])

        visualize_mols(frags, cols=2, rows=2, bond_idx=bond_des, colors=colors, wbos=wbos, fname='{}/{}_bond_{}_{}_frags_fixed.pdf'.format(name, name, bond_des[0], bond_des[1]),
                       align_to=2)




