"""
Visualize results of different fragmenting with different parameters. This script generates two figures for every
selected molecule in the benchmark set.
1. All fragments, with bond that fragment was built around highlighted with its ELF10 WBO. Titles ar parameters used
to get this fragment
2. Ridge plots to visualize WBO distributions for each fragment in results
"""

import fragmenter
from openeye import oechem, oedepict, oegraphsim, oequacpac

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sbn

import json


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


def to_pdf(molecules, bond_map_idx, fname, rows=3, cols=2, align=None):
    """
    Generate PDF of list of oemols or SMILES

    Parameters
    ----------
    molecules : list of OEMols
        These mols need to have map indices on bond of interest and WBO attached to that bond's data
    fname : str
        Name of PDF
    rows : int
        How many rows of molecules per page
    cols : int
        How many columns of molecule per page
    bond_map_idx : tuple of bond to highlight
    align: oemol
        molecule to align all other molecules in the list

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

    if align:
        if isinstance(align, str):
            ref_mol = oechem.OEGraphMol()
            oechem.OESmilesToMol(ref_mol, align)
        elif isinstance(align, (oechem.OEMol, oechem.OEMolBase, oechem.OEGraphMol)):
            ref_mol = align
        oedepict.OEPrepareDepiction(ref_mol)

    for i, mol in enumerate(molecules):
        cell = report.NewCell()
        mol_copy = oechem.OEMol(mol)
        oedepict.OEPrepareDepiction(mol_copy, False, True)

        atom_bond_set = oechem.OEAtomBondSet()
        a1 = mol_copy.GetAtom(oechem.OEHasMapIdx(bond_map_idx[0]))
        a2 = mol_copy.GetAtom(oechem.OEHasMapIdx(bond_map_idx[1]))
        b = mol_copy.GetBond(a1, a2)
        opts.SetBondPropertyFunctor(fragmenter.chemi.LabelWibergBondOrder())
        atom_bond_set.AddAtom(a1)
        atom_bond_set.AddAtom(a2)
        atom_bond_set.AddBond(b)
        hstyle = oedepict.OEHighlightStyle_BallAndStick
        hcolor = oechem.OEColor(oechem.OELightBlue)

        overlaps = oegraphsim.OEGetFPOverlap(ref_mol, mol_copy, oegraphsim.OEGetFPType(oegraphsim.OEFPType_Tree))
        oedepict.OEPrepareMultiAlignedDepiction(mol_copy, ref_mol, overlaps)
        disp = oedepict.OE2DMolDisplay(mol_copy, opts)
        oedepict.OEAddHighlighting(disp, hcolor, hstyle, atom_bond_set)

        oedepict.OERenderMolecule(cell, disp)
        oedepict.OEDrawCurvedBorder(cell, oedepict.OELightGreyPen, 10.0)

    oedepict.OEWriteReport(fname, report)


def ridge_plot(data, fname):
    """
    Generate Ridge plot of result WBO distributions
    Parameters
    ----------
    data: dict
        {bond_tuple: {parameters: {'frag': SMILES, 'wbo_dist': []}}}
    fname: str
        filename
    """
    with PdfPages(fname) as pdf:
        for bond in data:
            if bond == 'provenance':
                continue
            n = len(data[bond])
            fig = plt.figure()
            gs = gridspec.GridSpec(n, 2, width_ratios=[20, 1])
            fig.dpi = 400
            x_min = 3
            x_max = 0
            for p in data[bond]:
                wbos = data[bond][p]['wbo_dist']
                if min(wbos) < x_min:
                    x_min = min(wbos)
                if max(wbos) > x_max:
                    x_max = max(wbos)
            for i, p in enumerate(data[bond]):
                label = p.split('_')[0]
                wbos = data[bond][p]['wbo_dist']
                ax = plt.subplot(gs[i, :-1])
                ax.spines['left'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.spines['top'].set_visible(False)
                ax.spines['bottom'].set_visible(False)
                ax.patch.set_facecolor('none')
                if p == 'parent':
                    sbn.kdeplot(wbos, shade=True, color='red', alpha=1.0)
                    sbn.kdeplot(wbos, lw=3, color='black')
                else:
                    sbn.kdeplot(wbos, shade=True, color='steelblue', alpha=0.3)
                    sbn.kdeplot(wbos, lw=0.4, color='steelblue')
                if len(wbos) < 2:
                    for x in wbos:
                        plt.axvline(x=x, ymin=0, ymax=1, color='steelblue', linewidth=2.0, alpha=0.8)
                plt.xlim(x_min - 0.1, x_max + 0.1)
                plt.yticks([])
                ax.yaxis.set_label_coords(-0.05, 0)
                plt.ylabel(label, rotation=0, size=8)
                if i != n - 1:
                    plt.xticks([])
                else:
                    plt.xlabel('Wiberg Bond order')
                if i == 0:
                    plt.title(bond)
            # Magic to get overlapping distributions
            overlap = 0.5
            h_pad = 5 + (- 5 * (1 + overlap))
            fig.tight_layout(h_pad=h_pad)
            pdf.savefig(bbox_inches='tight')
            plt.close()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="calculate OE WBO")
    parser.add_argument('-n', '--name', type=str, help='Molecule name with number of its state')
    args = parser.parse_args()
    name = args.name

    with open('{}/{}_wbo_dists.json'.format(name, name), 'r') as f:
        results = json.load(f)
    with open('{}/{}_pfizer_wbo_dists.json'.format(name, name), 'r') as f:
        pfizer_results = json.load(f)

    mols = {}
    for bond in results:
        results[bond]['pfizer'] = pfizer_results[bond]
        if bond == 'provenance':
            continue
        print(bond)
        b = fragmenter.utils.deserialize_bond(bond)
        mols[b] = {'all': []}
        for parameters in results[bond]:

            mol = oechem.OEMol()
            oechem.OESmilesToMol(mol, results[bond][parameters]['frag'])
            mol.SetTitle(parameters)
            if not results[bond][parameters]['wbo_dist']:
                print('Not sure why wbo dist is missing. Recacluate it')
                charged = fragmenter.chemi.get_charges(mol, keep_confs=-1, strict_stereo=False, strict_types=False)
                for conf in charged.GetConfs():
                    mol_copy = oechem.OEMol(conf)
                    if oequacpac.OEAssignPartialCharges(mol_copy, oequacpac.OECharges_AM1BCCSym):
                        bo = get_bond(mol_copy, b)
                        wbo = bo.GetData('WibergBondOrder')
                        results[bond][parameters]['wbo_dist'].append(wbo)
            #charged_mol = fragmenter.chemi.get_charges(mol, strict_types=False, strict_stereo=False)
            #bo = get_bond(charged_mol, b)
            wbo = results[bond][parameters]['elf10_wbo']
            #results[bond][parameters]['elf10_wbo'] = wbo
            for a in mol.GetAtoms():
                if not a.GetMapIdx() in b:
                    a.SetMapIdx(0)
            bo = get_bond(mol, b)
            bo.SetData('WibergBondOrder', wbo)
            if parameters == 'parent':
                mols[b]['parent'] = mol
            mols[b]['all'].append(mol)

    #with open('{}/{}_wbos_dists_with_elf10.json'.format(name, name), 'w') as f:
    #    json.dump(results, f, indent=2, sort_keys=True)
    for b in mols:
        to_pdf(mols[b]['all'], bond_map_idx=b, align=mols[b]['parent'],
               fname='{}/{}_bond_{}_{}.pdf'.format(name, name, str(b[0]), str(b[1])))
    ridge_plot(results, fname='{}/{}_ridge.pdf'.format(name, name))

