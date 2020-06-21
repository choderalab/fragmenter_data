import json
import matplotlib.pyplot as plt
import seaborn as sbn
import numpy as np

from fragmenter.utils import deserialize_bond

sbn.set_style('whitegrid')
sbn.set_context('paper', font_scale=1.5)


def rmse(pre, tar):
    return np.sqrt((pre - tar)**2).mean()

with open('torsiondrive-job-idx.json', 'r') as f:
    job_indices = json.load(f)
with open('qcarchive_torsionscans.json', 'r') as f:
    td_drives = json.load(f)

angles = np.arange(-180, 195, 15)
for mol in job_indices:
    for bond in job_indices[mol]:
        plt.figure()
        b = deserialize_bond(bond)
        fname = '{}_bond_{}_{}.pdf'.format(mol, b[0], b[1])
        for i, frag in enumerate(['parent', 'pfizer', '0.03']):
            index = job_indices[mol][bond][frag]
            if index not in td_drives:
                continue
            if frag == 'pfizer':
                label = 'Simple'
            if frag == '0.03':
                label = 'WBO'
            if frag == 'parent':
                label = 'Parent'
            if frag in ['pfizer', '0.03']:
                target_index = job_indices[mol][bond]['parent']
                rms = rmse(np.asarray(td_drives[target_index]['energy']), np.asarray(td_drives[index]['energy']))
                label = label + ' RMSE: {}'.format(round(rms, 2))

            plt.plot(angles, td_drives[index]['energy'], '-', linewidth=1.0, color=sbn.color_palette('colorblind')[i])
            plt.plot(angles, td_drives[index]['energy'], '.', linewidth=1.0, color=sbn.color_palette('colorblind')[i], label=label)
        plt.legend()
        plt.xlabel('Torsion angle (degree)')
        plt.ylabel('Relative Energy (kJ/mol)')
        plt.title('QC torsion scans')
        plt.tight_layout()
        plt.savefig(fname)
        plt.close()


