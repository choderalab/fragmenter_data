import json
import glob
import matplotlib.pyplot as plt
import seaborn as sbn
import numpy as np

def mmd_x_xsqred(x, y):
    """
    Maximum mean discrepancy with squared kernel
    This will distinguish mean and variance
    see https://stats.stackexchange.com/questions/276497/maximum-mean-discrepancy-distance-distribution
    Parameters
    ----------
    x : list of ints
    y : list of ints

    Returns
    -------
    mmd score

    """

    y_arr = np.asarray(y)
    y_squared = y_arr*y_arr
    x_arr = np.asarray(x)
    x_squared = np.square(x_arr)

    E_x = np.mean(x_arr)
    E_y = np.mean(y_arr)

    E_x_squared = np.mean(x_squared)
    E_y_squared = np.mean(y_squared)

    mmd2 = (E_x - E_y)**2 + (E_x_squared - E_y_squared)**2
    return np.sqrt(mmd2)

names = glob.glob('*/')
scores_omega = []
scores_scans = []
scores_combined = []
differences = []
bonds = []
for name in names:
    name = name.split('/')[0]
    with open('{}/{}_wbo_dists_fixed.json'.format(name, name), 'r') as f:
        omega_results = json.load(f)
    with open('{}/{}_pfizer_wbo_dists.json'.format(name, name), 'r') as f:
        pfizer_results = json.load(f)
    with open('{}/{}_wbo_scans_fixed.json'.format(name, name), 'r') as f:
        scan_results = json.load(f)
    for bond in omega_results:
        if bond == 'provenance':
            continue
        bonds.append(name + '_' + bond)
        parent_omega_wbos = omega_results[bond]['parent']['wbo_dist']
        wbo_fragment_omega_wbos = omega_results[bond]['0.03']['wbo_dist']
        pfizer_fragment_omega_wbos = pfizer_results[bond]['wbo_dist']
        omega_score_wbo = mmd_x_xsqred(parent_omega_wbos, wbo_fragment_omega_wbos)
        omega_score_pfizer = mmd_x_xsqred(parent_omega_wbos, pfizer_fragment_omega_wbos)
        scores_omega.append((omega_score_pfizer, omega_score_wbo))

        parent_scan_wbos = scan_results[bond]['parent']['wbos']
        wbo_fragment_scan_wbos = scan_results[bond]['wbo_scheme']['wbos']
        pfizer_fragment_scan_wbos = scan_results[bond]['pfizer']['wbos']
        scan_score_wbo = mmd_x_xsqred(parent_scan_wbos, wbo_fragment_scan_wbos)
        scan_score_pfizer = mmd_x_xsqred(parent_scan_wbos, pfizer_fragment_scan_wbos)
        scores_scans.append((scan_score_pfizer, scan_score_wbo))

        parent_combined = parent_omega_wbos + parent_scan_wbos
        wbo_fragment_combined = wbo_fragment_omega_wbos + wbo_fragment_scan_wbos
        pfizer_fragment_combined = pfizer_fragment_omega_wbos + pfizer_fragment_scan_wbos
        combined_score_wbo = mmd_x_xsqred(parent_combined, wbo_fragment_combined)
        combined_score_pfizer = mmd_x_xsqred(parent_combined, pfizer_fragment_combined)
        scores_combined.append((combined_score_pfizer, combined_score_wbo))

combined_differences = [i-j for i,j in scores_combined]
plt.rcParams.update({'font.size': 12})
array = np.asarray(combined_differences)
x1 = array[(array <= 1) & (array > 0)]
x2 = array[~((array <= 1) & (array >= 0))]
x3 = array[(array == 0)]
counts, bins = np.histogram(array, bins=30)
# finally, do the plot
f, ax = plt.subplots()
sbn.distplot(x1, bins=bins, kde=False, rug=True, color=sbn.color_palette('colorblind')[0], label='WBO scheme better ({} fragments)'.format(len(x1)))
sbn.distplot(x2, bins=bins, kde=False, color=sbn.color_palette('colorblind')[4], rug=True, label='Pfizer scheme better ({} fragments)'.format(len(x2)))
sbn.distplot(x3, bins=bins, color=sbn.color_palette('colorblind')[2], label='Equally good ({} fragments)'.format(len(x3)))
plt.legend()
plt.xlim(-0.1);
plt.xlabel('Differences in distance score (Pfizer scheme - WBO scheme)')
plt.ylabel('Counts');
plt.savefig('combined-score-differences-fixed.pdf')
