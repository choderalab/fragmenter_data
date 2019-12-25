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
    ig = plt.figure(figsize=(8,8))
    gs = gridspec.GridSpec(3, 3)
    ax_main = plt.subplot(gs[1:3, :2])
    ax_xDist = plt.subplot(gs[0, :2])
    ax_yDist = plt.subplot(gs[1:3, 2])

    ax_main.grid(True)
    ax_main.scatter(x, y, alpha=0.7, edgecolor='black', zorder=2)
    ax_main.set(xlabel="Distance Score", ylabel="Computational cost")


    # Remove Nans
    xs = [i for i in x if not np.isnan(i)]
    kde = stats.gaussian_kde(xs)
    xx = np.linspace(-100, max(xs)+1, 100000)
    ax_xDist.plot(sorted(xx),kde(sorted(xx)), color='black')
    ax_xDist.set_yticks([])
    ax_xDist.tick_params(labelbottom=False)
    ax_xDist.set_xlim(-0.05, 0.4)
    ax_xDist.set_ylim(0, 80)
    ax_xDist.fill_betweenx(kde(sorted(xx)), 0, sorted(xx), alpha=0.3)
    ax_xDist.spines['left'].set_visible(False)
    ax_xDist.spines['right'].set_visible(False)
    ax_xDist.spines['top'].set_visible(False)
    ax_main.set_xlim(-0.05, 0.4)

    ys = [i for i in y if not np.isnan(i)]
    kde_y = stats.gaussian_kde(ys)
    yy = np.linspace(-100000, max(ys)+1, 100000)
    ax_yDist.plot(kde_y(sorted(yy)), sorted(yy), color='black')
    ax_yDist.fill_betweenx(sorted(yy), 0, kde_y(sorted(yy)), alpha=0.3)
    ax_yDist.set_xticks([])
    ax_yDist.tick_params(labelleft=False)
    ax_yDist.set_ylim(-5000, 70000)
    ax_yDist.set_xlim(0, 0.0003)
    ax_yDist.spines['top'].set_visible(False)
    ax_yDist.spines['right'].set_visible(False)
    ax_yDist.spines['bottom'].set_visible(False)
    ax_main.set_ylim(-5000, 70000)

    plt.subplots_adjust(wspace=0, hspace=0)
    plt.savefig(fname)


if __name__ == '__main__':
    scores = {}
    too_big = []
    too_bad = []
    names = glob.glob('*/')
    for n in names:
        n = n[:-1]
        with open('{}/{}_wbo_dists.json'.format(n, n), 'r') as f:
            wbos = json.load(f)
        for bond in wbos:
            if bond == 'provenance' or bond == 'p':
                continue
            for param in wbos[bond]:
                if param == 'parent':
                    continue
                params = param.split('_')
                if 'path' in params:
                    params.remove('length')
                threshold = float(params[0])
                hueristic = params[1]
                rotors = params[2]
                if rotors == 'True':
                    # Do not keep non rotor substituent unless it is meta to bond of interest. It just makes fragments bigger
                    # but does not improve their score.
                    continue
                # key = '{}_{}'.format(hueristic, rotors)
                if len(params) > 3:
                    f = params[-1]
                    if f == 'False':
                        # Do not use results from scheme that did not tag functional groups. It can have weird results and should not be done.
                        continue

                if threshold not in scores:
                    scores[threshold] = {}
                if hueristic not in scores[threshold]:
                    scores[threshold][hueristic] = {'scores': [], 'size': []}
                parent = wbos[bond]['parent']['wbo_dist']
                y = wbos[bond][param]['wbo_dist']
                score = mmd_x_xsqred(x=parent, y=y)
                heavy_atoms = n_heavy_atoms(wbos[bond][param]['frag'])
                if threshold == 0.05 and heavy_atoms > 25 and hueristic == 'path':
                    too_big.append((n, bond, wbos[bond][param]['frag']))
                if threshold == 0.01 and score > 0.2 and hueristic == 'path':
                    too_bad.append((n, bond, wbos[bond][param]['frag']))
                scores[threshold][hueristic]['scores'].append(score)
                scores[threshold][hueristic]['size'].append(heavy_atoms)
    print('Could not find a small enough fragment for:')
    print(too_big)
    print('Could not find a fragment with low enough score:')
    print(too_bad)

    # Plot distributions
    for i in (0.001, 0.005, 0.01, 0.05, 0.1):
        joint_plot(scores[i]['path']['scores'], np.asarray(scores[i]['path']['size']) ** 3,
                   fname='jointplot_{}.pdf'.format(i))


