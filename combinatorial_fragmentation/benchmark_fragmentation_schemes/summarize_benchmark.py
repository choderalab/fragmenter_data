"""
This script summarizes the combinatorial fragmentation, benchmarking experiment.
After finding the top 100 scoring bonds (bonds that generated fragments with very high distance scores. These bonds in
their corresponding molecules are challanging to fragment because of significant non-local effects. A few end up returning
the parent molecule because the effects are very long distance)

This script collects the results of different fragmentation scheme, calculated the distance scores (mmd with squared kernel)
and computational score (heavy_atoms**3 as an estimate O(n^3) of DFT) and generates joint plots for different thresholds
tested.

Note on other parameters
1. Not tagging functional groups listed in fragmenter data yaml file is a not a good idea. You can end up with weird fragments
2. There is no need to pull along non rotatable substituents. It just makes the fragments larger and no significant improevment on score
3. shortest path length is a better hueristic than greatest WBO. It leads to smaller fragments and similar scores.

There were some molecules where the scheme did not find the most optimal fragment. Some more analysis is needed for that
"""

from openeye import oechem
import json
import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sbn
import numpy as np
import glob
from scipy import stats


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

def n_heavy_atoms(smiles):
    """
    """
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, smiles)
    n = 0
    for a in mol.GetAtoms():
        if not a.IsHydrogen():
            n += 1
    return n


def joint_plot(x, y, fname):
    """
    A scatter plot with KDEs of marginals on the side
    Parameters
    ----------
    x : list
        values for x
    y : list
        values for y
    fname : str
        filename
    """

    #sbn.set_style('whitegrid')
    #sbn.set_context('paper', font_scale=1.7)
    plt.rcParams.update({'font.size': 14})
    ig = plt.figure(figsize=(8,8))
    gs = gridspec.GridSpec(3, 3)
    ax_main = plt.subplot(gs[1:3, :2])
    ax_xDist = plt.subplot(gs[0, :2])
    ax_yDist = plt.subplot(gs[1:3, 2])

    ax_main.grid(True)
    ax_main.scatter(x, y, alpha=0.5, edgecolor='black', zorder=2)
    ax_main.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8], minor=False)
    ax_main.set(xlabel="Distance Score", ylabel=r'CPU seconds $\propto$ NHeavy$^{2.6}$')


    # Remove Nans
    xs = [i for i in x if not np.isnan(i)]
    kde = stats.gaussian_kde(xs)
    xx = np.linspace(-100, max(xs)+1, 100000)
    ax_xDist.plot(sorted(xx),kde(sorted(xx)), color='black')
    ax_xDist.set_yticks([])
    ax_xDist.tick_params(labelbottom=False)
    ax_xDist.set_xlim(-0.05, 0.8, 0.1)
    ax_xDist.set_ylim(0, 80)
    ax_xDist.fill_betweenx(kde(sorted(xx)), 0, sorted(xx), alpha=0.3)
    ax_xDist.spines['left'].set_visible(False)
    ax_xDist.spines['right'].set_visible(False)
    ax_xDist.spines['top'].set_visible(False)
    ax_main.set_xlim(-0.05, 0.8, 0.1)

    ys = [i for i in y if not np.isnan(i)]
    kde_y = stats.gaussian_kde(ys)
    yy = np.linspace(-100000, max(ys)+1, 100000)
    ax_yDist.plot(kde_y(sorted(yy)), sorted(yy), color='black')
    ax_yDist.fill_betweenx(sorted(yy), 0, kde_y(sorted(yy)), alpha=0.3)
    ax_yDist.set_xticks([])
    ax_yDist.tick_params(labelleft=False)
    ax_yDist.set_ylim(-500, 15000)
    ax_yDist.set_xlim(0, 0.001)
    ax_yDist.spines['top'].set_visible(False)
    ax_yDist.spines['right'].set_visible(False)
    ax_yDist.spines['bottom'].set_visible(False)
    ax_main.set_ylim(-500, 15000)

    plt.subplots_adjust(wspace=0, hspace=0)
    plt.savefig(fname)


if __name__ == '__main__':
    scores = {}
    too_big = []
    too_bad = []
    pfizer_failure = []
    lower_left = {}
    names = glob.glob('*/')
    for n in names:
        n = n[:-1]
        with open('{}/{}_wbo_dists.json'.format(n, n), 'r') as f:
            wbos = json.load(f)
        with open('{}/{}_pfizer_wbo_dists.json'.format(n, n), 'r') as f:
            pfizer_results = json.load(f)
        for bond in wbos:
            wbos[bond]['pfizer'] = pfizer_results[bond]
            if bond == 'provenance' or bond == 'p':
                continue
            for threshold in wbos[bond]:
                # if param == 'parent':
                #     continue
                # params = param.split('_')
                # if not params[0] == 'pfizer':
                #     if 'path' not in params:
                #         continue
                #     else:
                #         params.remove('length')
                #     threshold = params[0]
                #     hueristic = params[1]
                #     rotors = params[2]
                #     if rotors == 'True':
                #         continue
                #     if len(params) > 3:
                #         f = params[3]
                #         if f == 'False':
                #             continue
                # if param == 'pfizer':
                #     threshold = 'pfizer'
                #     hueristic = 'pfizer'

                if threshold not in scores:
                    print(threshold)
                    scores[threshold] = {'scores':[], 'size': []}
                    lower_left[threshold] = {'lower_left': 0, 'outside': 0}
                # if hueristic not in scores[threshold]:
                #     print(hueristic)
                #     scores[threshold][hueristic] = {'scores': [], 'size': []}
                #     lower_left[threshold][hueristic] = {'lower_left': 0, 'outside': 0}
                parent = wbos[bond]['parent']['wbo_dist']
                y = wbos[bond][threshold]['wbo_dist']
                score = mmd_x_xsqred(x=parent, y=y)
                if 'frag' not in wbos[bond][threshold]:
                    print('frag not in dictionary')
                    print(n)
                    print(bond)
                    print(threshold)
                    print(wbos[bond][threshold].keys())
                heavy_atoms = n_heavy_atoms(wbos[bond][threshold]['frag'])

                if score < 0.05 and heavy_atoms**2.6 < 4000:
                    lower_left[threshold]['lower_left'] += 1
                else:
                    lower_left[threshold]['outside'] += 1
                if threshold == ('0.03', '0.05', '0.1') and heavy_atoms > 25 :
                    too_big.append((n, bond, wbos[bond][threshold]['frag']))
                if threshold in ('0.01', '0.03', '0.05') and score > 0.2 :
                    too_bad.append((n, bond, wbos[bond][threshold]['frag']))
                if threshold == 'pfizer' and score > 0.3:
                    pfizer_failure.append((n, bond, wbos[bond][threshold]['frag']))
                scores[threshold]['scores'].append(score)
                scores[threshold]['size'].append(heavy_atoms)
    print('Could not find a small enough fragment for:')
    print(too_big)
    print('Could not find a fragment with low enough score:')
    print(too_bad)
    print('pfizer failure')
    print(pfizer_failure)

    # Plot distributions
    print(scores.keys())
    for i in ('0.001', '0.005', '0.01', '0.03',  '0.05', '0.07',  '0.1'):
        joint_plot(scores[i]['scores'], np.asarray(scores[i]['size']) ** 2.6,
                   fname='jointplot_{}.pdf'.format(i))
    print(scores['pfizer'].keys())
    print(max(scores['pfizer']['scores']))
    joint_plot(scores['pfizer']['scores'], np.asarray(scores['pfizer']['size'])** 2.6,
               fname='jointplot_pfizer_fixed.pdf')
    with open('summary.json', 'w') as f:
        json.dump(lower_left, f, indent=2, sort_keys=True)


