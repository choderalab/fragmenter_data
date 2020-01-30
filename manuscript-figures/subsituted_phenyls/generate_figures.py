import json
import seaborn as sbn
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
import pandas as pd
import math

import qcfractal.interface as ptl
from fragmenter.utils import HARTREE_2_KJMOL
from fragmenter import chemi
from openeye import oedepict, oechem, oegraphsim

color_keys = ['maroon', 'brown', 'indianred', 'red', 'coral','orange', 'gold', 'darkkhaki','yellowgreen','limegreen',
              'mediumseagreen', 'teal', 'steelblue', 'cornflowerblue', 'royalblue', 'darkblue',
              'mediumblue', 'slateblue', 'blueviolet', 'purple','mediumvioletred', 'deeppink', 'hotpink',
              'palevioletred', 'pink', 'lightpink']

fgroup_symbols_colors = {
    #'phenoxide': 'C[O-]',
    'dimethylamino': (r'$\mathrm{\mathsf{N(Me)_2}}$', color_keys[0]),
    'methylamino': (r'$\mathrm{\mathsf{NHMe}}$', color_keys[1]),
    'amino': (r'$\mathrm{\mathsf{NH_2}}$', color_keys[2]),
    'ethylamino': (r'$\mathrm{\mathsf{NHEt}}$', color_keys[3]),
    'propylamino': (r'$\mathrm{\mathsf{NH(C_3H_7)}}$', color_keys[4]),
    'hydroxy': (r'$\mathrm{\mathsf{OH}}$', color_keys[5]),
    'methoxy': (r'$\mathrm{\mathsf{OMe}}$', color_keys[6]),
    'ethoxy': (r'$\mathrm{\mathsf{OEt}}$', color_keys[7]),
    'dimethylurea': (r'$\mathrm{\mathsf{NHCON(Me)_2}}$', color_keys[8]),
    'urea': (r'$\mathrm{\mathsf{NHCONHMe}}$', color_keys[9]),
    'phenylurea': (r'$\mathrm{\mathsf{NHCONH_2}}$', color_keys[10]),
    'ethylamide': (r'$\mathrm{\mathsf{NHCOEt}}$', color_keys[11]),
    'amide': (r'$\mathrm{\mathsf{NHCOMe}}$', color_keys[12]),
    'fluoro': (r'$\mathrm{\mathsf{F}}$', color_keys[13]),
    'chloro': (r'$\mathrm{\mathsf{Cl}}$', color_keys[14]),
    'cyano': (r'$\mathrm{\mathsf{CN}}$', color_keys[15]),
    'methyl': (r'$\mathrm{\mathsf{Me}}$', color_keys[16]),
    'bromo': (r'$\mathrm{\mathsf{Br}}$', color_keys[17]),
    'carbamate': (r'$\mathrm{\mathsf{OCONH_2}}$', color_keys[18]),
    'benzoicacid': (r'$\mathrm{\mathsf{COOH}}$', color_keys[19]),
    'iodo': (r'$\mathrm{\mathsf{I}}$', color_keys[20]),
    'ethoxycarbonyl': (r'$\mathrm{\mathsf{COOEt}}$', color_keys[21]),
    'trimethylamonium': (r'$\mathrm{\mathsf{N(Me)_3^+}}$', color_keys[22]),
    'trifluoromethyl': (r'$\mathrm{\mathsf{CF_3}}$', color_keys[23]),
    'nitro': (r'$\mathrm{\mathsf{NO_2}}$', color_keys[24])
}

# Generate joy plot
fgroup_wbos = {}
for fgroup in fgroup_symbols_colors:
    if fgroup not in fgroup_wbos:
        fgroup_wbos[fgroup] = []
    with open('../../phenyl_benchmark/data/{}_R1_wbos.json'.format(fgroup), 'r') as f:
        wbos = json.load(f)
    for w in wbos:
        fgroup_wbos[fgroup].append(w[0])

colors = mcolors.CSS4_COLORS

fig, axes = plt.subplots(len(fgroup_wbos))
for i, fgroup in enumerate(fgroup_wbos):
    ax = plt.subplot(len(fgroup_wbos), 1, i+1)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.patch.set_facecolor('none')
    sbn.kdeplot(fgroup_wbos[fgroup], shade=True, alpha=0.6,
                color=colors[fgroup_symbols_colors[fgroup][1]])
    sbn.kdeplot(fgroup_wbos[fgroup], shade=False, color='black', lw=0.8)
    plt.xlim(0.70, 1.4)
    plt.yticks([])
    ax.yaxis.set_label_coords(-0.05, 0)
    plt.ylabel(fgroup_symbols_colors[fgroup][0], rotation=0, size=8,
               color=colors[fgroup_symbols_colors[fgroup][1]])
    if i == len(fgroup_wbos)-1:
        plt.xlabel('Bond order')
    else:
        plt.xticks([])

overlap=1.0
h_pad = 5 + (- 5*(1 + overlap))
fig.tight_layout(h_pad=h_pad)
plt.savefig('figures/wbo_dist_joy_plot.pdf')


# See if there is a correlation with Hammet sigma parameters. Values were taken from
# doi:10.1021/cr00002a004
subs = ['H','dimethylamino', 'methylamino', 'amino', 'ethylamino', 'hydroxy', 'methoxy', 'phenylurea', 'amide',
        'fluoro', 'chloro','cyano', 'methyl', 'bromo', 'benzoicacid', 'ethoxycarbonyl', 'trifluoromethyl', 'nitro']
sigma_m = [0.0, -0.16, -0.21, -0.16, -0.24, 0.12, 0.12, -0.02, 0.21, 0.34, 0.37, 0.56, -0.07, 0.39, 0.37, 0.37, 0.43, 0.71]
sigma_p = [0.0, -0.83, -0.70, -0.66, -0.61, -0.37, -0.27, -0.24, 0.0, 0.06, 0.23, 0.66, -0.17, 0.45, 0.45, 0.45,  0.54, 0.78]
wbo_cooh_meta = [0.96, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.96, 0.96, 0.95, 0.95, 0.95, 0.96, 0.95, 0.96, 0.95, 0.95, 95]
wbo_cooh_para = [0.96, 0.97, 0.97, 0.97, 0.97, 0.96, 0.96, 0.97, 0.97, 0.96, 0.96, 0.96, 0.96, 0.96, 0.95, 0.95, 0.95, 95]
wbo_r_meta = [0.96, 1.07, 1.08, 1.12, 1.08, 1.06, 1.04, 1.02, 1.02, 1.02, 1.0,  1.0, 1.0, 0.99, 0.96, 0.93, 0.91, 0.85]
wbo_r_para = [0.96, 1.11, 1.10, 1.12, 1.14, 1.08, 1.05, 1.04, 1.03, 1.03, 1.01, 1.0, 1.0, 0.99, 0.95, 0.93, 0.91, 0.85]

hammet_sigmas = {'substituent':subs, 'sigma_p': sigma_p, 'sigma_m': sigma_m, 'wbo_cooh_meta': wbo_cooh_meta,
                 'wbo_cooh_para': wbo_cooh_para,'wbo_r_meta': wbo_r_meta, 'wbo_r_para': wbo_r_para}
df = pd.DataFrame(hammet_sigmas)

# plot correlation
markersize=9
fontsize=6
for sigma in ('m', 'p'):
    fig, ax = plt.subplots()
    for row in df.iterrows():
        if sigma == 'm':
            x = row[1].wbo_r_meta
            y = row[1].sigma_m
        if sigma == 'p':
            x = row[1].wbo_r_para
            y = row[1].sigma_p
        if row[1].substituent == 'H':
            plt.plot(x, y, '.', color='black', markersize=markersize, label='H')
            plt.annotate('H', (x, y),
                     textcoords='offset points', xytext=(3, 2), color= 'black', fontsize=fontsize)
            continue
        plt.plot(x, y, '.', markersize=markersize, color=fgroup_symbols_colors[row[1].substituent][1],
                     label=fgroup_symbols_colors[row[1].substituent][0])
        plt.annotate(fgroup_symbols_colors[row[1].substituent][0], (x, y),
                     textcoords='offset points', xytext=(3, 2), color= fgroup_symbols_colors[row[1].substituent][1], fontsize=fontsize)

    plt.xlim(0.83, 1.16)
    plt.ylim(-0.86, 0.85)
    plt.ylabel(r'$\sigma_{}$'.format(sigma))
    plt.xlabel('Wiberg Bond Order');
    if sigma == 'm':
        r_value = df.corr().sigma_m.wbo_r_meta
    if sigma == 'p':
        r_value = df.corr().sigma_p.wbo_r_para
    print(r_value)
    textstr = r'$\rho =%.2f$' % (r_value)
    props = dict(boxstyle='square', facecolor='white', alpha=0.5)
    ax.text(0.75, 0.95, textstr, transform=ax.transAxes, fontsize=12,
            verticalalignment='top', bbox=props)
    plt.tight_layout()
    fig.savefig('figures/hammett_sigma_{}.pdf'.format(sigma))


# Generate torsion barrier height vs ELF10 AM1 WBO plot
with open('../../phenyl_benchmark/data/qcarchive_torsiondrives.json', 'r') as f:
    fgroups_td = json.load(f)

fig, ax = plt.subplots()
for fgroup in fgroup_symbols_colors:
    if fgroup not in fgroups_td:
        continue
    energies = fgroups_td[fgroup]['energy']
    am1_wbos = fgroups_td[fgroup]['elf10_am1_wbo']
    max_energies = [max(energy) for energy in energies]
    slope, intercept, r_value, p_value, std_err = stats.linregress(am1_wbos, max_energies)
    fgroups_td[fgroup]['stats'] = [slope, r_value**2, p_value, std_err]
    plt.plot(np.unique(am1_wbos), np.poly1d([slope, intercept])(np.unique(am1_wbos)), fgroup_symbols_colors[fgroup][1], label=fgroup_symbols_colors[fgroup][0])
    plt.scatter(x=am1_wbos, y=max_energies, color=fgroup_symbols_colors[fgroup][1], s=4)

l = ax.legend(bbox_to_anchor=(1, 1))
for line, text in zip(l.get_lines(), l.get_texts()):
    text.set_color(line.get_color())

plt.xlabel('ELF10 AM1 Wiberg bond order')
plt.ylabel('Torsion barrier height (kJ/mol)')
plt.tight_layout()
plt.savefig('figures/energy_vs_wbo.pdf')

# generate table
stats_table = {'functional group': [], 'slope': [], 'r^2': [], 'P value': [], 'standard error': []}
for fgroup in fgroup_symbols_colors:
    if fgroup not in fgroups_td:
        continue
    stats_table['functional group'].append(fgroup_symbols_colors[fgroup][0])
    stats_table['slope'].append(fgroups_td[fgroup]['stats'][0])
    stats_table['r^2'].append(fgroups_td[fgroup]['stats'][1])
    stats_table['P value'].append(fgroups_td[fgroup]['stats'][2])
    stats_table['standard error'].append(fgroups_td[fgroup]['stats'][3])
latex_table = pd.DataFrame(stats_table).to_latex(index=False)
with open('figures/stats_ov.tex', 'w') as f:
    f.write(latex_table)