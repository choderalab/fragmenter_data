import json
from openeye import oechem, oequacpac
import glob

def oemol_from_old_json(tagged_smiles, xyz_string):
    """
    Generate oemol from old json input files. Since I am having trouble generating the same conformations as
    generated for the first hf3c study, and I want to do a 1 to 1 comparison to see if AM1 can reproduce
    hf3c correlations, I am using the old input files. With several verison updates with openeye and openmoltools,
    it's hard to get the same conformers.
    """
    oemol = oechem.OEMol()
    oechem.OESmilesToMol(oemol, tagged_smiles)

    mapping = {a.GetIdx(): a.GetMapIdx() for a in oemol.GetAtoms()}
    symbols = []
    xyz = []
    lines = xyz_string.split('\n')
    for l in lines:
        line = l.split()
        if len(line) <=1:
            continue
        symbols.append(line[0])
        xyz.append([float(i) for i in line[1:]])
    xyz_mapped = []
    for i in range(len(symbols)):
        xyz_mapped.extend(xyz[mapping[i]-1])

    oemol.SetCoords(oechem.OEFloatArray(xyz_mapped))
    oemol.SetDimension(3)

    return oemol


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Calculate AM1 WBOs')
    parser.add_argument('-n', '--name', type=str,
                        help='name of kinase inhibitor')

    args = parser.parse_args()
    name = args.name
    files = glob.glob('../data_generation/data/{}/*input.json'.format(name))

    oemols = []
    for file in files:
        with open(file, 'r') as f:
            im = json.load(f)
        oemol = oemol_from_old_json(im['tagged_smiles'], im['molecule'])
        oemols.append(oemol)

    # Calculate AM1 WBO and save oemols (this takes ~ .5 hour)
    ofs = oechem.oemolostream()
    ofs.open('{}_for_am1_comparison.oeb'.format(name))
    oemols_wbo = []
    for i, mol in enumerate(oemols):
        print(i)
        mol_copy = oechem.OEMol(mol)
        if oequacpac.OEAssignPartialCharges(mol_copy, oequacpac.OECharges_AM1BCCSym):
            oemols_wbo.append(mol_copy)
            oechem.OEWriteMolecule(ofs, mol_copy)
            pass
        else:
            print('am1 failed for {}'.format(i))

    print('calculated am1 wbos for {} molecules'.format(len(oemols_wbo)))
