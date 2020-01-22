"""
Download torsiondrive data for figure
Ran on 2020-1-21.
QCPortal version 0.13.0
openeye version 2019.Oct.2
"""
import qcportal as ptl
import json
from openeye import oechem

client = ptl.FractalClient()
collections = ['OpenFF Group1 Torsions', 'SMIRNOFF Coverage Torsion Set 1']
for c in collections:

    td_dataset = client.get_collection('TorsionDriveDataset', c)
    print(td_dataset.status(['default']))
    opts_per_td = {'heavy_atoms': [], 'opts_per_td': []}
    gradients_per_opts = {'heavy_atoms': [], 'gradients_per_opt': []}
    dictionary = {'driver': [],'method': [], 'basis': [], 'nbasis': [], 'nalpha': [], 'nbeta': [],
                  'natoms': [], 'heavy_atoms': [], 'cpu': [], 'nthreads': [], 'wall_time': [] }

    for i, index in enumerate(td_dataset.df.index):
        print(i, index)
        entry = td_dataset.get_entry(index)
        smiles = entry.attributes['canonical_isomeric_smiles']
        mol = oechem.OEMol()
        oechem.OESmilesToMol(mol, smiles)
        heavy_atoms = sum([1 for i in mol.GetAtoms() if i.GetAtomicNum() > 1])
        td = td_dataset.get_record(index, 'default')
        history = td.get_history()
        n_opts = 0
        for angle in history:
            for opt in history[angle]:
                n_opts += 1
                opt_traj = opt.get_trajectory()
                gradients_per_opts['gradients_per_opt'].append(len(opt_traj))
                gradients_per_opts['heavy_atoms'].append(heavy_atoms)
                for t in opt_traj:
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
        opts_per_td['opts_per_td'].append(n_opts)
        opts_per_td['heavy_atoms'].append(heavy_atoms)

        with open('off_torsions_cpu_time.json', 'w') as f:
            json.dump(dictionary, f, indent=2, sort_keys=True)
        with open('opts_per_td.json', 'w') as f:
            json.dump(opts_per_td, f, indent=2, sort_keys=True)
        with open('gradients_per_opt.json', 'w') as f:
            json.dump(gradients_per_opts, f, indent=2, sort_keys=True)


