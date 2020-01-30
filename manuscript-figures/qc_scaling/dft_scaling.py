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

drug_mols = pd.read_csv('../drugbank_small_mols.csv')

# Select on CPU to be consistent

df_intel = df[df['cpu'].str.contains('Intel*')]
df_amd = df[df['cpu'].str.contains('AMD*')]
df_intel_one = df[df['cpu'] == 'Intel(R) Xeon(R) CPU E5-2697 v4 @ 2.30GHz']

# Fit to a power curve to get order of scale
def power_law(x, m, c):
    return  x**m * c

def compute_stat(x, y):
    return curve_fit(power_law, x, y, method='dogbox')[0]

df_intel_one['cpu_time'] = df_intel_one['wall_time']*df_intel_one['nthreads']
x = df_intel_one['heavy_atoms'].values
y = df_intel_one['cpu_time'].values

popt, pcov = curve_fit(power_law, x, y, method='dogbox')
# Boostrap 95% CI
# def compute_stat(x, y):
#     return curve_fit(power_law, x, y, method='dogbox')[0]
ci = arch.bootstrap.IIDBootstrap(x, y).conf_int(compute_stat, 1000)

# Now do a bunch of manipulations to get everything to line up right
# This is done to trick seaborn to leave spaces where we are missing data
heavy_atoms = [0, 1, 2, 3, 4, 20, 21, 22, 24, 25, 30, 32, 33, 34, 35, 36]
wall_time = [-100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100]

sbn.set_style('whitegrid')
sbn.set_context('paper', font_scale=1.5)

fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(8, 10))

# Fitted power law
ax[0].plot(list(sorted(x)), power_law(sorted(x), *popt), color='black')
ax[0].plot(list(sorted(x)), power_law(sorted(x), *ci[0]), '--', color='grey')
ax[0].plot(list(sorted(x)), power_law(sorted(x), *ci[1]), '--', color='grey')
ax[0].fill_between(sorted(x), power_law(sorted(x), *ci[0]), power_law(sorted(x), *ci[1]), alpha=0.2, color='grey')
textstr = r'CPU time = %.2f*(NHeavy)$^{%.2f}$' % (popt[1], popt[0])
props = dict(boxstyle='square', facecolor='white', alpha=0.5)
ax[0].text(0.05, 0.95, textstr, transform=ax[0].transAxes, fontsize=14,
        verticalalignment='top', bbox=props);

# Boxplot
unique = np.sort(df_intel_one["heavy_atoms"].unique())
x_new = list(df_intel_one['heavy_atoms'].values) + heavy_atoms
y_new = list(df_intel_one['wall_time']*df_intel_one['nthreads']) + wall_time
sbn.boxplot(x_new, y_new, fliersize=0.3, ax=ax[0])
ax[0].set_ylim(0, 4200)
xlables = [' ', ' ', ' ', ' ', ' ', '5', ' ', '7', ' ', '9', ' ', '11', ' ', '13', ' ', '15', ' ', '17', ' ', '19', ' ',
          ' ', ' ', '23', ' ', ' ', '26', ' ', '28', ' ', ' ', '31', ' ', ' ', ' ', ' ', ' ', '37', ' ']
ax[0].set_xticks(unique, minor=True);
ax[0].xaxis.grid(True, which='minor')
ax[0].xaxis.set_ticks_position('bottom')
ax[0].tick_params(which='major', width=0.75, length=2.5)
ax[0].tick_params(which='minor', width=1.0, length=5.0)
ax[0].set_xlim(0, 38)
ax[0].set_xticklabels([])
ax[0].set_ylabel('CPU time (seconds)')


# Heavy atoms distribution
sbn.kdeplot(drug_mols.heavy_atoms, shade=True, color='grey', ax=ax[1], legend=False)
ax[1].set_xlim(0, 38)
xlables = ['5', ' ', '7', ' ', '9', ' ', '11', ' ', '13', ' ', '15', ' ', '17', ' ', '19', '23', '26', ' ', '28', ' ', '31', '37']
ax[1].set_xticks(unique, minor=False);
ax[1].xaxis.grid(True, which='minor')
ax[1].xaxis.set_ticks_position('bottom')
ax[1].tick_params(which='major', width=0.75, length=2.5)
ax[1].tick_params(which='minor', width=1.0, length=5.0)
ax[1].set_xticklabels(xlables)
#ax[1].set_yticklabels([])
ax[1].set_ylabel('Fraction of DrugBank')
ax[1].set_xlabel('Heavy atoms')

ax[0].set_title('CPU time as a funtion of heavy atoms');
ax[1].set_title('Distribution of DrugBank small molecule sizes')
plt.savefig('B3LYP_scaling.pdf', bbox_inches='tight')

# Generate SI figure. Separate for Intel and AMD CPUs
def find_missing(lst):
    return [x for x in range(0, 38)
            if x not in lst]

sbn.set_context('paper', font_scale=1.3)
cpu_seen = []
fig, ax = plt.subplots(nrows=4, ncols=3, figsize=(12, 16))
i = 0
j = 0
n = 0
for cpu in df_intel['cpu']:
    if cpu == 'Intel(R) Xeon(R) CPU E5-2697 v4 @ 2.30GHz':
        # this one is in main text. No need for SI
        continue
    if cpu in cpu_seen:
        continue
    if not 'Xeon(R)' in cpu:
        print(cpu)
        print('This should not happen')
    # print(cpu)
    n += 1
    print(i, j)
    new_df = df[(df['cpu'] == cpu)]
    unique = np.sort(new_df["heavy_atoms"].unique())
    missing = find_missing(unique)
    fake_data = []
    for value in missing:
        fake_data.append(-1000)

    x_new = list(new_df['heavy_atoms']) + missing
    y_new = list(new_df['wall_time'] * new_df['nthreads']) + fake_data

    xlabels = ['', '1', '', '', '', '5', '', '', '', '9', '', '', '', '13', '', '', '', '17', '',
               '', '', '21', '', '', '', '25', '', '', '', '29', '', '', '', '33', '', '', '', '37']
    sbn.boxplot(x_new, y_new, fliersize=0.5, ax=ax[i, j])
    ax[i, j].set_xticks(unique, minor=True);
    ax[i, j].xaxis.grid(True, which='minor')
    ax[i, j].xaxis.set_ticks_position('bottom')
    ax[i, j].tick_params(which='major', width=0.75, length=2.5)
    ax[i, j].tick_params(which='minor', width=1.0, length=5.0)
    if i == 3:
        ax[i, j].set_xticklabels(xlabels)
        ax[i, j].set_xlabel('Heavy atoms');
    else:
        ax[i, j].set_xticklabels([])

    ax[i, j].set_ylim(0, 4500)
    ax[i, j].set_xlim(0, 38)
    if len(unique) > 2:
        new_df['cpu_time'] = new_df['wall_time'] * new_df['nthreads']
        x = new_df['heavy_atoms'].values
        y = new_df['cpu_time'].values

        popt, pcov = curve_fit(power_law, x, y, method='dogbox')
        ci = arch.bootstrap.IIDBootstrap(x, y).conf_int(compute_stat, 1000)
        ax[i, j].plot(list(sorted(x_new)), power_law(sorted(x_new), *popt), color='black')
        textstr = r'%.2f*(NHeavy)$^{%.2f}$' % (popt[1], popt[0])
        props = dict(boxstyle='square', facecolor='white', alpha=0.0)
        ax[i, j].text(0.02, 0.98, textstr, transform=ax[i, j].transAxes,  # fontsize=14,
                      verticalalignment='top', bbox=props);
        ax[i, j].plot(list(sorted(x_new)), power_law(sorted(x_new), *ci[0]), '--', color='grey')
        ax[i, j].plot(list(sorted(x_new)), power_law(sorted(x_new), *ci[1]), '--', color='grey')
        ax[i, j].fill_between(sorted(x_new), power_law(sorted(x_new), *ci[0]), power_law(sorted(x_new), *ci[1]),
                              alpha=0.2, color='grey')

    if j == 0:
        ax[i, j].set_ylabel('CPU time (seconds)')
    else:
        ax[i, j].set_yticklabels([])

    ax[i, j].set_title(cpu[17:]);
    cpu_seen.append(cpu)
    j += 1
    if (n % 3) == 0:
        j = 0
        i += 1
plt.savefig('SI_Intel_scaling.pdf', bbox_inches='tight')

# AMD
sbn.set_context('paper', font_scale=1.5)
cpu_seen = []
fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(10, 10))
print(ax)
i = 0
j = 0
n = 0
for cpu in df_amd['cpu']:
    if cpu == 'Intel(R) Xeon(R) CPU E5-2697 v4 @ 2.30GHz':
        # this one is in main text. No need for SI
        continue
    if cpu in cpu_seen:
        continue
    n += 1
    print(i, j)
    new_df = df[(df['cpu'] == cpu)]
    unique = np.sort(new_df["heavy_atoms"].unique())
    missing = find_missing(unique)
    fake_data = []
    for value in missing:
        fake_data.append(-10000)

    x_new = list(new_df['heavy_atoms']) + missing
    y_new = list(new_df['wall_time'] * new_df['nthreads']) + fake_data

    xlabels = ['', '1', '', '', '', '5', '', '', '', '9', '', '', '', '13', '', '', '', '17', '',
               '', '', '21', '', '', '', '25', '', '', '', '29', '', '', '', '33', '', '', '', '37']
    sbn.boxplot(x_new, y_new, fliersize=0.5, ax=ax[i, j])
    ax[i, j].set_xticks(unique, minor=True);
    ax[i, j].xaxis.grid(True, which='minor')
    ax[i, j].xaxis.set_ticks_position('bottom')
    ax[i, j].tick_params(which='major', width=0.75, length=2.5)
    ax[i, j].tick_params(which='minor', width=1.0, length=5.0)
    # ax[i, j].set_xticks([1, 5, 9, 13, 17, 21, 25, 29, 33, 37], minor=False);
    if i == 1:
        ax[i, j].set_xticklabels(xlabels)
        ax[i, j].set_xlabel('Heavy atoms');
    else:
        ax[i, j].set_xticklabels([])

    ax[i, j].set_ylim(0, 20000)
    ax[i, j].set_xlim(0, 38)

    if j == 0:
        ax[i, j].set_ylabel('CPU time (seconds)')
    else:
        ax[i, j].set_yticklabels([])

    ax[i, j].set_title(cpu);
    cpu_seen.append(cpu)
    j += 1
    if (n % 2) == 0:
        j = 0
        i += 1
plt.savefig('SI_AMD_scaling.pdf', bbox_inches='tight')

# Plot gradients per optimization and optimizations per torsiondrive
dfs = []
with open('gradients_per_opt_ki.json', 'r') as f:
    dfs.append(pd.read_json(f))
with open('gradients_per_opt.json', 'r') as f:
    dfs.append(pd.read_json(f))
with open('gradients_per_opt_0.json', 'r') as f:
    dfs.append(pd.read_json(f))

for i in range(1, 13):
    with open('gradients_per_opt_0_{}.json'.format(i), 'r') as f:
        dfs.append(pd.read_json(f))

df = pd.concat(dfs)

fig, ax = plt.subplots()
unique = np.sort(df["heavy_atoms"].unique())

missing = find_missing(unique)
fake_data = []
for value in missing:
    fake_data.append(-1000)

x = df['heavy_atoms']
y = df['gradients_per_opt']
#popt, pcov = curve_fit(power_law, x, y, method='dogbox')

x_new = list(df['heavy_atoms']) + missing
y_new = list(df['gradients_per_opt']) + fake_data
# ax.plot(list(sorted(x_new)), power_law(sorted(x_new), *popt), color='black')
# ci = arch.bootstrap.IIDBootstrap(x, y).conf_int(compute_stat, 100)
# textstr = r'%.2f*(NHeavy)$^{%.2f}$' % (popt[1], popt[0])
# props = dict(boxstyle='square', facecolor='white', alpha=0.80)
# ax.text(0.53, 0.98, textstr, transform=ax.transAxes,  # fontsize=14,
#         verticalalignment='top', bbox=props);
# ax.plot(list(sorted(x_new)), power_law(sorted(x_new), *ci[0]), '--', color='grey')
# ax.plot(list(sorted(x_new)), power_law(sorted(x_new), *ci[1]), '--', color='grey')
# ax.fill_between(sorted(x_new), power_law(sorted(x_new), *ci[0]), power_law(sorted(x_new), *ci[1]), alpha=0.2,
#                 color='grey')

xlabels = ['', '1', '', '', '', '5', '', '', '', '9', '', '', '', '13', '', '', '', '17', '',
           '', '', '21', '', '', '', '25', '', '', '', '29', '', '', '', '33', '', '', '', '37']
sbn.boxplot(x=x_new, y=y_new, fliersize=0.3, ax=ax)
ax.set_ylim(0)
ax.set_xticks(unique, minor=True);
ax.xaxis.grid(True, which='minor')
ax.xaxis.set_ticks_position('bottom')
ax.tick_params(which='major', width=0.75, length=2.5)
ax.tick_params(which='minor', width=1.0, length=5.0)
ax.set_xticklabels(xlabels);
ax.set_xlabel('Heavy atoms')
ax.set_ylabel('Number of gradient evaluations')
plt.title('Gradient evaluations per optimization');
plt.savefig('gradients_per_opts.pdf', bbox_inches='tight')

dfs = []
with open('opts_per_td.json', 'r') as f:
    dfs.append(pd.read_json(f))
with open('opts_per_td_0.json', 'r') as f:
    dfs.append(pd.read_json(f))

for i in range(1, 13):
    with open('opts_per_td_0_{}.json'.format(i), 'r') as f:
        dfs.append(pd.read_json(f))
df = pd.concat(dfs)

def find_missing(lst):
    return [x for x in range(0, 20)
            if x not in lst]

fig, ax = plt.subplots()
unique = np.sort(df["heavy_atoms"].unique())

missing = find_missing(unique)
fake_data = []
for value in missing:
    fake_data.append(-1000)

x = df['heavy_atoms']
y = df['opts_per_td']
popt, pcov = curve_fit(power_law, x, y, method='dogbox')

x_new = list(df['heavy_atoms']) + missing
y_new = list(df['opts_per_td']) + fake_data
# ax.plot(list(sorted(x_new)), power_law(sorted(x_new), *popt), color='black')
# ci = arch.bootstrap.IIDBootstrap(x, y).conf_int(compute_stat, 100)
# textstr = r'%.2f*(NHeavy)$^{%.2f}$' % (popt[1], popt[0])
# props = dict(boxstyle='square', facecolor='white', alpha=0.80)
# ax.text(0.53, 0.98, textstr, transform=ax.transAxes, #fontsize=14,
# verticalalignment='top', bbox=props);
# ax.plot(list(sorted(x_new)), power_law(sorted(x_new), *ci[0]), '--', color='grey')
# ax.plot(list(sorted(x_new)), power_law(sorted(x_new), *ci[1]), '--', color='grey')
# ax.fill_between(sorted(x_new), power_law(sorted(x_new), *ci[0]), power_law(sorted(x_new), *ci[1]), alpha=0.2, color='grey')

xlabels = ['', '1', '', '', '', '5', '', '', '', '9', '', '', '', '13', '', '', '', '17', '',
           '', '', '21', '', '', '', '25', '', '', '', '29', '', '', '', '33', '', '', '', '37']
sbn.boxplot(x=x_new, y=y_new, fliersize=0.3, ax=ax)
ax.set_ylim(0)
ax.set_xticks(unique, minor=True);
ax.xaxis.grid(True, which='minor')
ax.xaxis.set_ticks_position('bottom')
ax.tick_params(which='major', width=0.75, length=2.5)
ax.tick_params(which='minor', width=1.0, length=5.0)
ax.set_xticklabels(xlabels);
plt.title('Optimizations per torsion drive')
ax.set_xlabel('Heavy atoms')
ax.set_ylabel('Number of optimizations per torsion drive')
plt.savefig('opts_per_td.pdf', bbox_inches='tight')