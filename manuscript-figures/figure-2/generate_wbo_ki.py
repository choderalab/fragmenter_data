import fragmenter
import pandas as pd
from openeye import oechem
import time


kinase_inhibitors = fragmenter.chemi.file_to_oemols('kinase_inhibitors.smi')

initial_time = time.time()

nmolecules = 0
ofs = oechem.oemolostream()
if not ofs.open('kinase_inhibitors_wbo.oeb'):
    oechem.OEThrow.Fatal('Unable to open file for writing')
for i, sm in enumerate(kinase_inhibitors):
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, sm)
    try:
        charged_mol = fragmenter.chemi.get_charges(mol, strict_stereo=False, strict_types=False)
    except RuntimeError:
        print('charging failed for {}'.format(sm))
    try:
        oechem.OEWriteMolecule(ofs, charged_mol)
    except Exception as e:
        print('Failed to charge {}'.format(sm))
        print(e)
    nmolecules += 1
    total_time = time.time() - initial_time
    average_time = total_time / nmolecules
    print(f'{nmolecules} molecules processed in {total_time} seconds : {average_time} seconds/molecule')
