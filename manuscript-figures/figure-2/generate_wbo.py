import fragmenter
import pandas as pd
from openeye import oechem
import time
import sys

def main(argv=[__name__]):
    print(argv)
    if len(argv) != 4:
        oechem.OEThrow.Usage('"%s <infile> <outfile> <frag> <nfrags>" % argv[0]')

    frag = int(argv[1])
    nfrags = int(argv[2])
    small_mols = pd.read_csv('drugbank_small_mols.csv')

    filtered_drugbank = small_mols.loc[(small_mols['heavy_atoms'] <= 50) &
                                       (small_mols['fda_approved'] == True) &
                                       (small_mols['connected_components'] == 1)]


    # open file for writing
    nmolecules = len(filtered_drugbank.smiles)
    nstart = int( (nmolecules / nfrags) * (frag-1) )
    nprocess = min( nmolecules, int( (nmolecules / nfrags) * frag ) - int( (nmolecules / nfrags) * (frag-1) ) )
    print(f'Fragment {frag} of {nfrags} : Starting at molecule {nstart} and processing {nprocess} molecules to write to {argv[3]}')
    initial_time = time.time()

    nmolecules = 0
    ofs = oechem.oemolostream()
    if not ofs.open(argv[3]):
        oechem.OEThrow.Fatal('Unable to open file for writing')
    for i, sm in enumerate(filtered_drugbank.smiles):
        if (i < nstart) or (i >= nstart+nprocess):
            continue
        mol = oechem.OEMol()
        oechem.OESmilesToMol(mol, sm)
        try:
            charged_mol = fragmenter.chemi.get_charges(mol, strict_stereo=False, strict_types=False)
        except RuntimeError:
            print('charging failed for {}'.format(sm))
        oechem.OEWriteMolecule(ofs, charged_mol)
        nmolecules += 1
        total_time = time.time() - initial_time
        average_time = total_time / nmolecules
        print(f'{nmolecules} molecules processed in {total_time} seconds : {average_time} seconds/molecule')

if __name__ == "__main__":
    sys.exit(main(sys.argv))