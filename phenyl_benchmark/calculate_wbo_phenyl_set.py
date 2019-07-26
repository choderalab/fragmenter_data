"""
Calculate WBO for R1 on phenyl ring as R2 changes and generate joy plot
"""

import matplotlib.pyplot as plt
from openeye import oechem
from fragmenter import chemi
import seaborn as sbn
import json
import matplotlib.colors as mcolors


def group_by_fgroup_and_wbo(fgroup, molecules):
    """
    Group molecules by functional group at the R1 position and calculate the WBO
    for R1 in different chemical environments
    Parameters
    ----------
    fgroup : str
        functional group
    molecules : oemol

    Returns
    -------
    wbo_dict: dictionary of functional groups to WBO distribution at R1

    """
    mols = []
    for mol in molecules:
        name = mol.GetTitle().split('_')
        if fgroup in name:
            mols.append(mol)

    # remove duplicates
    smiles = []
    deduplicated_mols = []
    for mol in mols:
        sm = oechem.OEMolToSmiles(mol)
        if sm not in smiles:
            smiles.append(sm)
            deduplicated_mols.append(mol)
    # sort by wbo
    wbo_dict = {}
    for mol in deduplicated_mols:
        name = mol.GetTitle().split('_')
        charged = chemi.get_charges(mol)
        if fgroup == name[-1]:
            map_idx = [1]
        elif fgroup == name[1]:
            map_idx = [2, 3]
        for bond in charged.GetBonds():
            m1 = bond.GetBgn().GetMapIdx()
            m2 = bond.GetEnd().GetMapIdx()
            if (m1 in map_idx or m2 in map_idx) and not bond.IsInRing():
                wbo = bond.GetData('WibergBondOrder')
                if wbo not in wbo_dict:
                    wbo_dict[wbo] = []
                wbo_dict[wbo].append(charged)
    return wbo_dict

def prep_mols_for_vis(wbo_dict, fgroup):
    # first sort wbo
    wbo_keys = sorted(list(wbo_dict.keys()))
    molecules = []
    bond_map_idx = [(4,5)]*len(wbo_keys)
    wbos = []
    for bo in wbo_keys:
        drop = False
        mol = oechem.OEMol(wbo_dict[bo][0])
        name = mol.GetTitle().split('_')
        if fgroup == name[-1]:
            map_idx = [1]
        elif fgroup == name[1]:
            map_idx = [2, 3]
        for bond in mol.GetBonds():
            a1 = bond.GetBgn()
            a2 = bond.GetEnd()
            m1 = a1.GetMapIdx()
            m2 = a2.GetMapIdx()
            if name[-1] == fgroup:
                to_check_ortho = [1, 4, 5]
            elif name[1] == fgroup:
                to_check_ortho = [2, 3, 4, 5]
            if (m1 in to_check_ortho or m2 in to_check_ortho) and bond.IsInRing():
                if a1.GetAtomicNum() == 7 or a2.GetAtomicNum() == 7:
                    # The N is ortho to the bond. drop out from list
                    drop = True
            if (m1 in map_idx or m2 in map_idx) and not bond.IsInRing():
                bond.GetBgn().SetMapIdx(4)
                bond.GetEnd().SetMapIdx(5)
                wbo = bond.GetData('WibergBondOrder')
            else:
                # Remove wbo so that only R1 WBO is generated in visualization
                to_delete = bond.GetData('WibergBondOrder')
                tag = oechem.OEGetTag('WibergBondOrder')
                bond.DeleteData(tag)
        if not drop:
            molecules.append(mol)
            wbos.append(wbo)
    return molecules, bond_map_idx, wbos


# Load all molecules
all_mols = chemi.file_to_oemols('phenyls.smi')
all_mols.extend(chemi.file_to_oemols('pyridine_ortho.smi'))
all_mols.extend(chemi.file_to_oemols('pyridine_meta.smi'))
all_mols.extend(chemi.file_to_oemols('pyridine_para.smi'))

"""
Only look at R1 functional groups that have torsions (skip phenoxice, fluoro, chloro, cyano, bromo and iodo).
This leaves 20 functional groups.
"""
fgroups_smarts = {
    #'phenoxide': 'C[O-]',
    'dimethylamino': 'CN(C)(C)',
    'methylamino': 'CNC',
    'amino': 'CN',
    'ethylamino': 'CNCC',
    'propylamino': 'CNCCC',
    'hydroxy': 'CO',
    'methoxy': 'COC',
    'ethoxy': 'COCC',
    'dimethylurea': 'CNC(=O)N(C)(C)',
    'urea': 'CNC(=O)NC',
    'phenylurea': 'CNC(=O)N',
    'ethylamide': 'CNC(=O)CC',
    'amide': 'CNC(=O)C',
    #'fluoro': 'CF',
    #'chloro': 'CCl',
    'methyl': 'CC',
    #'cyano': 'CC#N',
    #'bromo': 'CBr',
    'carbamate': 'COC(=O)N',
    #'iodo': 'CI',
    'benzoicacid': 'C(=O)O',
    'ethoxycarbonyl': 'CC(=O)OCC',
    'trifluoromethyl': 'CC(F)(F)(F)',
    'trimethylamonium': 'C[N+](C)(C)C',
    'nitro': 'C[N+](=O)[O-]'
}

color_keys = ['rosybrown', 'indianred', 'red', 'orange', 'gold', 'yellow','greenyellow', 'green', 'limegreen',
          'lightseagreen', 'teal', 'cyan', 'deepskyblue', 'royalblue', 'mediumslateblue', 'blueviolet', 'mediumorchid', 'violet',
          'palevioletred', 'lightpink']
colors = mcolors.CSS4_COLORS

#color=cm.rainbow_r(np.linspace(0,1,len(fgroups_smarts)))
#color = chemi._KELLYS_COLORS
def get_rgb_int(rgb, alpha):
    return (int(rgb[0]*255), int(rgb[1]*255), int(rgb[2]*255), int(rgb[3]*alpha))

fgroups_wbos = {}
for i, fgroup in enumerate(fgroups_smarts):
    print(fgroup)
    wbo_dict = group_by_fgroup_and_wbo(fgroup, all_mols)
    fgroups_wbos[fgroup] = wbo_dict
    mols, bond_maps, wbos = prep_mols_for_vis(wbo_dict, fgroup=fgroup)
    #c = get_rgb_int(color[i], alpha=180)
    c_oe = oechem.OEColor(colors[color_keys[i]])
    chemi.to_pdf(molecules=mols, bond_map_idx=bond_maps, bo=wbos, align=mols[0], color=c_oe,
    fname='figures/oe_wbos/{}.pdf'.format(fgroup))

    # Save WBOs
    to_save = [[wbo, mol.GetTitle(), oechem.OEMolToSmiles(mol)] for wbo, mol in zip(wbos, mols)]
    with open('data/{}_R1_wbos.json'.format(fgroup), 'w') as f:
        json.dump(to_save, f, indent=2, sort_keys=True)

# Generate joy plot
fig, axes = plt.subplots(len(fgroups_wbos))
for i, fgroup in enumerate(fgroups_wbos):
    ax = plt.subplot(len(fgroups_wbos), 1, i+1)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.patch.set_facecolor('none')
    sbn.kdeplot(list(fgroups_wbos[fgroup].keys()), shade=True, alpha=0.8, color=colors[color_keys[i]])
    sbn.kdeplot(list(fgroups_wbos[fgroup].keys()), shade=False, color='black', lw=1.0)
    plt.xlim(0.70, 1.4)
    plt.yticks([])
    ax.yaxis.set_label_coords(-0.05, 0)
    plt.ylabel(fgroup, rotation=0, size=8)
    if i == len(fgroups_wbos)-1:
        plt.xlabel('Bond order')
    else:
        plt.xticks([])

overlap=1.0
h_pad = 5 + (- 5*(1 + overlap))
fig.tight_layout(h_pad=h_pad)
plt.savefig('figures/oe_wbos/phenyl_set_wbo_dist.pdf')
