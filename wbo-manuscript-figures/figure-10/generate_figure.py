import json
from openeye import oechem, oedepict
import glob

mol_names = glob.glob('../../combinatorial_fragmentation/rank_fragments/selected/*')

pdf_filename = 'selected.pdf'

itf = oechem.OEInterface()
PageByPage = True
suppress_h = True
rows = 5
cols = 5
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

pen = oedepict.OEPen(oechem.OEBlack, oechem.OEBlack, oedepict.OEFill_Off, 1.0)
opts.SetDefaultBondPen(pen)
oedepict.OESetup2DMolDisplayOptions(opts, itf)

mols = []
for i, name in enumerate(mol_names):
    n = name.split('/')[-1]
    with open('{}/{}_selected_bonds.json'.format(name, n), 'r') as f:
        selected = json.load(f)
    smiles = selected['parent_smiles']
    bonds = selected['bonds']
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, smiles)
    atom_bond_set = oechem.OEAtomBondSet()
    for bond in bonds:
        a1 = mol.GetAtom(oechem.OEHasMapIdx(int(bond[0])))
        a2 = mol.GetAtom(oechem.OEHasMapIdx(int(bond[1])))
        atom_bond_set.AddAtom(a1)
        atom_bond_set.AddAtom(a2)
        b = mol.GetBond(a1, a2)
        atom_bond_set.AddBond(b)
    hstyle = oedepict.OEHighlightStyle_BallAndStick
    hcolor = oechem.OEColor(oechem.OEYellowTint)

    mol.SetTitle('')
    cell = report.NewCell()
    #mol_copy = oechem.OEMol(mol)
    oedepict.OEPrepareDepiction(mol, False, suppress_h)
    mols.append((mol, atom_bond_set))

# minscale = float("inf")
# for mol in mols:
#     minscale = min(minscale, oedepict.OEGetMoleculeScale(mol[0], opts))
#
# opts.SetScale(minscale)

for idx, cell in enumerate(report.GetCells()):
    disp = oedepict.OE2DMolDisplay(mols[idx][0], opts)
    oedepict.OEAddHighlighting(disp, hcolor, hstyle, mols[idx][1])
    oedepict.OERenderMolecule(cell, disp)
    #oedepict.OEDrawCurvedBorder(cell, oedepict.OELightGreyPen, 10.0)

oedepict.OEWriteReport(pdf_filename, report)
