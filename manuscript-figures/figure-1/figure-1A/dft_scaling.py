"""
script to generate scaling figure. Data used to generate this was the first 689 molecules in OpenFF Group1
Torsions TorsionDriveDataset (that was all I managed to download)
and the Kinase Inhibitors: WBO Distributions OptimizationDataset.

"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sbn
from scipy.optimize import curve_fit
import numpy as np
import arch.bootstrap


with open('ki_opt_cpu_time.json', 'r') as f:
    ki_stats = pd.read_json(f)
with open('off_torsions_cpu_time.json', 'r') as f:
    off_td_stats = pd.read_json((f))

dfs = [ki_stats, off_td_stats]
for i in range(1, 13):
       with open('off_torsions_cpu_time_0_{}.json'.format(i), 'r') as f:
            dfs.append(pd.read_json(f))
df = pd.concat(dfs)

# Select on CPU to be consistent

df_intel = df[df['cpu'].str.contains('Intel*')]
df_amd = df[df['cpu'].str.contains('AMD*')]
df_intel_one = df[df['cpu'] == 'Intel(R) Xeon(R) CPU E5-2697 v4 @ 2.30GHz']

# Fit to a power curve to get order of scale
def power_law(x, m, c):
    return  x**m * c

def compute_stat(x, y):
    return curve_fit(power_law, x, y, method='dogbox')[0][0]

df_intel_one['cpu_time'] = df_intel_one['wall_time']*df_intel_one['nthreads']
x = df_intel_one['heavy_atoms'].values
y = df_intel_one['cpu_time'].values

popt, pcov = curve_fit(power_law, x, y, method='dogbox')
# Boostrap 95% CI
def compute_stat(x, y):
    return curve_fit(power_law, x, y, method='dogbox')[0]
ci = arch.bootstrap.IIDBootstrap(x, y).conf_int(compute_stat, 1000)

# Now do a bunch of manipulations to get everything to line up right
# This is done to trick seaborn to leave spaces where we are missing data
heavy_atoms = [0, 1, 2, 3, 4, 20, 21, 22, 24, 25, 30, 32, 33, 34, 35, 36]
wall_time = [-100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100]

sbn.set_style('whitegrid')
sbn.set_context('paper', font_scale=1.5)

fig, ax = plt.subplots()
unique = np.sort(df_intel_one["heavy_atoms"].unique())
x_new = list(df_intel_one['heavy_atoms'].values) + heavy_atoms
y_new = list(df_intel_one['wall_time']*df_intel_one['nthreads']) + wall_time
sbn.boxplot(x_new, y_new, fliersize=0.3)
ax.set_ylim(0, 4500)
xlables = [' ', ' ', ' ', ' ', ' ', '5', ' ', '7', ' ', '9', ' ', '11', ' ', '13', ' ', '15', ' ', '17', ' ', '19', ' ',
          ' ', ' ', '23', ' ', ' ', '26', ' ', '28', ' ', ' ', '31', ' ', ' ', ' ', ' ', ' ', '37', ' ']
ax.set_xticks(unique, minor=True);
ax.xaxis.grid(True, which='minor')
ax.set_xticklabels(xlables)
ax.set_xlabel('Heavy atoms')
ax.set_ylabel('CPU time (seconds)')

plt.plot(list(sorted(x)), power_law(sorted(x), *popt), color='black')
plt.plot(list(sorted(x)), power_law(sorted(x), *ci[0]), '--', color='grey')
plt.plot(list(sorted(x)), power_law(sorted(x), *ci[1]), '--', color='grey')
plt.fill_between(sorted(x), power_law(sorted(x), *ci[0]), power_law(sorted(x), *ci[1]), alpha=0.2, color='grey')
textstr = r'CPU time = %.2f*(NHeavy)$^{%.2f}$' % (popt[1], popt[0])
props = dict(boxstyle='square', facecolor='white', alpha=0.5)
ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14,
        verticalalignment='top', bbox=props);
ax.set_title('CPU time as a funtion of heavy atoms');
plt.savefig('B3LYP_scaling.pdf', bbox_inches='tight')