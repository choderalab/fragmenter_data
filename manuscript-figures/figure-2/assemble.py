#!/usr/bin/env python3

"""
Stitch oeb files together
"""

import os
from openeye import oechem

# Determine prefixes to assembl
import glob
filenames = glob.glob("*.oeb")


print(f'Assembling dataset and deleting original fragments...')
filename = 'drugbank_am1_wbo.oeb'
ofs = oechem.oemolostream()
if not ofs.open(filename):
    oechem.OEThrow.Fatal("Unable to open %s for writing" % filename)

nmolecules = 0
for filename in filenames:

    if not os.path.exists(filename):
        continue

    ifs = oechem.oemolistream()
    if not ifs.open(filename):
        oechem.OEThrow.Fatal("Unable to open %s for reading" % filename)

    for mol in ifs.GetOEMols():
        oechem.OEWriteMolecule(ofs, mol)
        nmolecules += 1

    ifs.close()
    os.unlink(filename)
print(f'{nmolecules} molecules')

ofs.close()