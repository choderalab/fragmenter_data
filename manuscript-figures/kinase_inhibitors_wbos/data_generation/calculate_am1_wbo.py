import fragmenter
import argparse
import cmiles
from openeye import oequacpac, oechem
import time


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run psi4 to calculate bond orders')
    parser.add_argument('-n', '--name', type=str,
                        help='kinase inhibitor')

    args = parser.parse_args()
    name = args.name

    # load conformers
    conformers = fragmenter.chemi.file_to_oemols('data/{}.mol2'.format(name))

    ofs = oechem.oemolostream()
    ofs.open('data/{}_am1_wbo.oeb'.format(name))
    time_1 = time.time()
    for mol in conformers:
        cmiles.utils.add_atom_map(mol, in_place=True)
        if oequacpac.OEAssignPartialCharges(mol, oequacpac.OECharges_AM1BCCSym):
            oechem.OEWriteMolecule(ofs, mol)
        else:
            print('AM1 failed')
    time_2 = time.time()
    print('Calculated WBO for {} molecules in {} seconds'.format(len(conformers), time_2-time_1))

