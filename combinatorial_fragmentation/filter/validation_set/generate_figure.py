import pandas as pd
from openeye import oechem, oedepict


df = pd.read_csv('drugbank_filtered.csv')
smiles = df['smiles']
# Generate PDF of molecules in set
pdf_filename = 'filtered_db_figure.pdf'

itf = oechem.OEInterface()
PageByPage = True
suppress_h = True
rows = 10
cols = 10
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

pen = oedepict.OEPen(oechem.OEBlack, oechem.OEBlack, oedepict.OEFill_Off, 0.6)
opts.SetDefaultBondPen(pen)
oedepict.OESetup2DMolDisplayOptions(opts, itf)
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
opts.SetScale(minscale*2)

for mol in mols[::2]:
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, s)
    cell = report.NewCell()
    oedepict.OEPrepareDepiction(mol, False, suppress_h)
    #opts.SetTitleLocation(oedepict.OETitleLocation_Hidden)

for idx, cell in enumerate(report.GetCells()):
    #i = (len(mols)-1) - idx
    disp = oedepict.OE2DMolDisplay(mols[idx], opts)
    oedepict.OERenderMolecule(cell, disp)

oedepict.OEWriteReport(pdf_filename, report)
