"""
Download optimization data for figure
Ran on 2020-1-21.
QCPortal version 0.13.0
openeye version 2019.Oct.2
"""
import json
import qcportal as ptl
from openeye import oechem

client = ptl.FractalClient()
ki_ds = client.get_collection('OptimizationDataset', 'Kinase Inhibitors: WBO Distributions')

dictionary = {'driver': [],'method': [], 'basis': [], 'nbasis': [], 'nalpha': [], 'nbeta': [],
              'natoms': [], 'heavy_atoms': [], 'cpu': [], 'nthreads': [], 'wall_time': [] }

gradients_per_opts = {'heavy_atoms': [], 'gradients_per_opt': []}

for j, index in enumerate(ki_ds.df.index):
    print(j)
    entry = ki_ds.get_entry(index)
    smiles = entry.attributes['canonical_isomeric_smiles']
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, smiles)
    heavy_atoms = sum([1 for i in mol.GetAtoms() if i.GetAtomicNum() > 1])

    try:
        rec = ki_ds.get_record(index, 'default')
    except KeyError:
        continue
    traj = rec.get_trajectory()
    gradients_per_opts['heavy_atoms'].append(heavy_atoms)
    gradients_per_opts['gradients_per_opt'].append(len(traj))
    for t in traj:
        dictionary['driver'].append(t.driver.name)
        dictionary['method'].append(t.method)
        dictionary['basis'].append(t.basis)
        dictionary['nbasis'].append(t.properties.calcinfo_nmo)
        dictionary['natoms'].append(t.properties.calcinfo_natom)
        dictionary['nalpha'].append(t.properties.calcinfo_nalpha)
        dictionary['nbeta'].append(t.properties.calcinfo_nbeta)
        dictionary['cpu'].append(t.provenance.cpu)
        dictionary['nthreads'].append(t.provenance.nthreads)
        dictionary['wall_time'].append(t.provenance.wall_time)
        dictionary['heavy_atoms'].append(heavy_atoms)

with open('ki_opt_cpu_time.json', 'w') as f:
    json.dump(dictionary, f, indent=2, sort_keys=True)
with open('gradients_per_opt_ki.json', 'w') as f:
    json.dump(gradients_per_opts, f, indent=2, sort_keys=True)