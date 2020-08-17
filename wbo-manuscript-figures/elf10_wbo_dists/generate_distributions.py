"""
Generate ELF10 WBO distiributions for set of drug molecules.
---------
Ran on 2019-01-21
fragmenter version 0.0.4+24.g18bdc15.dirty
"""
import fragmenter
from openeye import oechem, oedepict
import matplotlib.pyplot as plt
import seaborn as sbn

mols = fragmenter.chemi.file_to_oemols('druglike_wbos.oeb')

# Sort bond by types to generate individual distributions
wbos = {'C-C': [],
        'ring': [],
        'H': [],
        'N': [],
        'O': [],
        'halognes': [],
        'P': [],
        'S': [],
        'all': [],
        'all_no_h': []}
missing_wbos = set()
for m in mols:
    for bond in m.GetBonds():
        if 'WibergBondOrder' not in bond.GetData():
            missing_wbos.add(m)
            continue
        wbo = bond.GetData('WibergBondOrder')
        wbos['all'].append(wbo)
        a1 = bond.GetBgn()
        a2 = bond.GetEnd()
        sym1 = oechem.OEGetAtomicSymbol(a1.GetAtomicNum())
        sym2 = oechem.OEGetAtomicSymbol(a2.GetAtomicNum())

        if bond.IsInRing():
            wbos['ring'].append(bond.GetData('WibergBondOrder'))
        if sym1 == 'N' or sym2 == 'N':
            wbos['N'].append(wbo)
        if sym1 == 'S' or sym2 == 'S':
            wbos['S'].append(wbo)
        if sym1 == 'P' or sym2 == 'P':
            wbos['P'].append(wbo)
        if sym1 == 'O' or sym2 == 'O':
            wbos['O'].append(wbo)
        if (sym1 == 'C' and sym2 == 'C') and not bond.IsInRing() and not sym1 == 'H' and not sym2 == 'H':
            wbos['C-C'].append(wbo)
        if sym1 == 'H' or sym2 == 'H':
            wbos['H'].append(wbo)
        if not sym1 == 'H' and not sym2 == 'H':
            wbos['all_no_h'].append(wbo)

plt.figure()
sbn.kdeplot(wbos['all'], shade=True, color=sbn.color_palette('colorblind')[0])
plt.xlabel('Wiberg Bond Order', fontsize=14)
plt.yticks([])
plt.xlim(0.1)
plt.xticks(fontsize=14)
plt.savefig('wbo_dist_all.pdf', bbox_inches='tight')

plt.figure()
sbn.kdeplot(wbos['C-C'],shade=True, label='C-C bonds (not in ring)', color=sbn.color_palette('colorblind')[0])
sbn.kdeplot(wbos['ring'], shade=True, label='bonds in rings', color=sbn.color_palette('colorblind')[4])
plt.xlabel('Wiberg Bond Order', fontsize=14)
plt.xticks(fontsize=14)
plt.yticks([])
plt.xlim(0.1)
plt.legend(fontsize=12)
plt.savefig('wbo_dist_carbon.pdf', bbox_inches='tight')

plt.figure()
hfont = {'fontname':'Helvetica'}
sbn.kdeplot(wbos['N'],shade=True, label='bonds with nitrogen', color=sbn.color_palette('colorblind')[0])
sbn.kdeplot(wbos['O'], shade=True, label='bonds with oxygen', color=sbn.color_palette('colorblind')[4])
plt.xlabel('Wiberg Bond Order', fontsize=14)
plt.yticks([])
plt.xticks(fontsize=14)
plt.legend(fontsize=12)
plt.xlim(0.1)
plt.savefig('wbo_dist_n_o.pdf', bbox_inches='tight')

# Generate PDF of molecules in set
pdf_filename = 'druglike_mols_for_wbo_dist.pdf'

itf = oechem.OEInterface()
PageByPage = True
suppress_h = True
rows = 9
cols = 7
ropts = oedepict.OEReportOptions(rows, cols)
ropts.SetHeaderHeight(0.001)
ropts.SetFooterHeight(0.001)
ropts.SetCellGap(0.001)
ropts.SetPageMargins(0.001)
report = oedepict.OEReport(ropts)
cellwidth, cellheight = report.GetCellWidth(), report.GetCellHeight()
opts = oedepict.OE2DMolDisplayOptions(cellwidth, cellheight, oedepict.OEScale_AutoScale)
opts.SetAromaticStyle(oedepict.OEAromaticStyle_Circle)
#opts.SetAtomColorStyle(oedepict.OEAtomColorStyle_WhiteMonochrome)

pen = oedepict.OEPen(oechem.OEBlack, oechem.OEBlack, oedepict.OEFill_Off, 0.5)
opts.SetDefaultBondPen(pen)
oedepict.OESetup2DMolDisplayOptions(opts, itf)
for m in mols[::-1]:
    cell = report.NewCell()
    oedepict.OEPrepareDepiction(m, False, suppress_h)
    opts.SetTitleLocation(oedepict.OETitleLocation_Hidden)

for idx, cell in enumerate(report.GetCells()):
    i = (len(mols)-1) - idx
    disp = oedepict.OE2DMolDisplay(mols[i], opts)
    oedepict.OERenderMolecule(cell, disp)

oedepict.OEWriteReport(pdf_filename, report)
