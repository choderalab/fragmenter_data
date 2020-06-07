import qcfractal.interface as ptl
import json

from fragmenter.utils import HARTREE_2_KJMOL
import numpy as np

def get_energies_wbo(dataset, index):
    """
    Collect energies and Lowdin-Wiberg bond orders from torsiondrive dataset

    Parameters
    ----------
    dataset: QCArchive TorsionDrive dataset
    index: str
        Job index

    Returns
    -------
    final_energies: list of final energy
    wbos: list of WBOs
    """
    angles = np.arange(-165, 195, 15)
    td = dataset.get_record(index, specification='default')
    final_energies =[td.get_final_energies(int(i))*HARTREE_2_KJMOL for i in angles]
    final_energies.insert(0, final_energies[-1])
    final_energies = list(np.asarray(final_energies) - min(final_energies))

    wbos = []
    n_atoms = len(td.get_final_molecules(-30).symbols)
    dih = td.keywords.dihedrals[0]
    for a in angles:
        opt = td.get_history(int(a), minimum=True)
        result = opt.get_trajectory()[-1]
        wiberg = np.array(result.extras['qcvars']['WIBERG_LOWDIN_INDICES']).reshape(-1, n_atoms)
        wbos.append(wiberg[dih[1], dih[2]])
    return final_energies, wbos

def get_conformers(dataset, index):
    """
    Download qcschema mols of optimized geometry
    """
    angles = np.arange(-165, 195, 15)
    td = dataset.get_record(index, specification='default')
    final_molecules =[td.get_final_molecules(int(i)).dict(encoding='json') for i in angles]
    return final_molecules

client = ptl.FractalClient()
validation_dataset = client.get_collection('TorsionDriveDataset', 'OpenFF Fragmenter Validation 1.0')

td_scans = {}

for index in validation_dataset.df.index:
    print(index)
    if 'I' in index:
        print('Job has iodine: {}'.format(index))
        continue
    try:
        final_energy, wbos = get_energies_wbo(validation_dataset, index)
    except KeyError:
        print('{} not complete'.format(index))
        continue
    td_scans[index] = {}
    td_scans[index]['energy'] = final_energy
    td_scans[index]['wbo'] = wbos

with open('qcarchive_torsionscans.json', 'w') as f:
    json.dump(td_scans, f, indent=2, sort_keys=True)
