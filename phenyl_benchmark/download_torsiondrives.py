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
phenyl_dataset = client.get_collection('TorsionDriveDataset', 'OpenFF Substituted Phenyl Set 1')

fgroups =  [
    'dimethylamino',
    'methylamino',
    'amino',
    'ethylamino',
    'propylamino',
    'hydroxy',
    'methoxy',
    'ethoxy',
    'dimethylurea',
    'urea',
    'phenylurea',
    'ethylamide',
    'amide',
    'carbamate',
    'benzoicacid',
    'ethoxycarbonyl',
    'nitro']

fgroups_td_scans = {}
fgroups_td_conformers = {}
for fgroup in fgroups:
    print(fgroup)
    fgroups_td_scans[fgroup] = {'indices': [], 'elf10_am1_wbo': [], 'energy': [], 'lowdin_wbos': []}
    fgroups_td_conformers[fgroup] = {}
    with open('data/{}_td_job_indices.json'.format(fgroup), 'r') as f:
        indices = json.load(f)
        for i in indices:
            if phenyl_dataset.get_record(i[0], specification='default').status == 'COMPLETE':
                if i[0] in fgroups_td_scans[fgroup]['indices']:
                    continue
                fgroups_td_scans[fgroup]['indices'].append(i[0])
                fgroups_td_scans[fgroup]['elf10_am1_wbo'].append(i[3])
                if i[0] in fgroups_td_conformers[fgroup]:
                    continue
                fgroups_td_conformers[fgroup][i[0]] = get_conformers(phenyl_dataset, i[0])
        with open('data/{}_qcarchive_td_conformers.json'.format(fgroup), 'w') as f:
            json.dump(fgroups_td_conformers, f, indent=2, sort_keys=True)

for fgroup in fgroups:
    print(fgroup)
    for index in fgroups_td_scans[fgroup]['indices']:
        energies, wbos = get_energies_wbo(phenyl_dataset, index)
        fgroups_td_scans[fgroup]['energy'].append(energies)
        fgroups_td_scans[fgroup]['lowdin_wbos'].append(wbos)

with open('data/qcarchive_torsiondrives.json', 'w') as f:
    json.dump(fgroups_td_scans, f, indent=2, sort_keys=True)