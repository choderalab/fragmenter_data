import matplotlib.pyplot as plt
import seaborn
import json
import glob
import numpy as np

# Get rescore JSON because that includes all fragmnents - not only the ones that have all 1-5 bonds. We need to see
# all kinds of chemical environment changes, even the ones that are not remote
files = glob.glob('../../combinatorial_fragmentation/rank_fragments/selected/*/rescore/*_oe_wbo_with_score.json')

dists_over_confs = []
dists_over_frags = []
for f in files:
    with open(f, 'r') as j:
        bond_wbos = json.load(j)
    for bond in bond_wbos:
        elf_dist = []
        for frag in bond_wbos[bond]:
            elf_estimate = bond_wbos[bond][frag]['elf_estimate']
            conf_dep_wbos = bond_wbos[bond][frag]['individual_confs']
            dists_over_confs.append(conf_dep_wbos)
            elf_dist.append(elf_estimate)
        dists_over_frags.append(elf_dist)

# Calculate stds for all distributions
stds_confs = [np.std(np.asarray(i)) for i in dists_over_confs]
stds_frags = [np.std(np.asarray(i)) for i in dists_over_frags]

seaborn.distplot(stds_confs, kde=True, color=seaborn.color_palette('colorblind')[0], label='conformations')
seaborn.distplot(stds_frags, kde=True, color=seaborn.color_palette('colorblind')[4], label='chemical environments')
plt.xlim(-0.002)
plt.legend(fontsize=14)
plt.yticks([])
plt.xticks(fontsize=14)
plt.xlabel('Standard deviation', fontsize=14)
plt.savefig('starndard_deviations.pdf', bbox_inches='tight')