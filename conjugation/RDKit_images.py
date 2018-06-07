
# coding: utf-8

# In[1]:


from IPython.display import SVG
from rdkit.Chem import AllChem as Chem 
from rdkit.Chem.Draw import rdMolDraw2D, MolToFile #IPythonConsole #Needed to show molecules
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions
import pandas as pd
from cairosvg import svg2png

from rdkit.Chem.Lipinski import RotatableBondSmarts
# In[29]:


# Generate images for all kinase inhibitors for both OE conjugation and RDKit conjugation
rdkit_mols = Chem.SmilesMolSupplier('kinase_inhibitors.smi')
for mol in rdkit_mols:
    # Write out image of resonance structure
    Chem.Compute2DCoords(mol)
    conjugated_bonds = [bond.GetIdx() for bond in mol.GetBonds() if bond.GetIsConjugated()]
    drawer = rdMolDraw2D.MolDraw2DSVG(800,400)
    drawer.DrawMolecule(mol, highlightAtoms=[], highlightBonds=conjugated_bonds)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText().replace('svg:','')
    outfile = open("images/{}_rdkit_res.png".format(mol.GetProp("_Name")), 'wb')
    svg2png(svg, write_to=outfile)
    outfile.close()
    
    # Find rotatable bonds
    atom_tuples = mol.GetSubstructMatches(RotatableBondSmarts)
    rotatable_bonds = []
    for a1, a2 in atom_tuples:
        bond = mol.GetBondBetweenAtoms(a1, a2)
        rotatable_bonds.append(bond.GetIdx())
    drawer = rdMolDraw2D.MolDraw2DSVG(800,400)
    drawer.DrawMolecule(mol, highlightAtoms=[], highlightBonds=rotatable_bonds)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText().replace('svg:','')
    outfile = open("images/{}_rdkit_rotor.png".format(mol.GetProp("_Name")), 'wb')
    svg2png(svg, write_to=outfile)
    outfile.close()

