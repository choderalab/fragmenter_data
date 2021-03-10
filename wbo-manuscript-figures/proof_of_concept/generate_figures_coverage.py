import json
import seaborn as sbn
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
import pandas as pd
import arch.bootstrap
import math
import qcfractal.interface as ptl
from fragmenter.utils import HARTREE_2_KJMOL
from fragmenter import chemi
from openeye import oedepict, oechem, oegraphsim
from openforcefield.topology import Molecule, Topology
from openforcefield.typing.engines.smirnoff import ForceField
import pickle

def checkTorsion(smiles, torsion_indices, ff_name):
    """
    Take mollist and check if the molecules in a list match a specific torsion id

        Parameters
        ----------
        molList : List of objects
            List of oemols with datatags generated in genData function

        Returns
        -------
        molList : list of objects
            List of oemol objects that have a datatag "IDMatch" that contain the torsion id
            involved in the QCA torsion drive
    """

    matches = []
    count = 0
    mols = []
    #tid=''
    #molecule = Molecule.from_mapped_smiles(smiles)
    print(smiles)
    from openeye import oechem
    # create a new molecule
    #mol = oechem.OEGraphMol()
    # convert the SMILES string into a molecule
    #oechem.OESmilesToMol(mol,smiles)
    #molecule = Molecule.from_smiles(smiles)
    #molecule=Molecule.from_openeye(mol)

    molecule = Molecule.from_mapped_smiles(smiles)
    topology = Topology.from_molecules(molecule)
    # Let's label using the Parsley force field
    forcefield = ForceField(ff_name, allow_cosmetic_attributes=True)
    # Run the molecule labeling
    molecule_force_list = forcefield.label_molecules(topology)
    params = []
    indices=[]
    # Print out a formatted description of the torsion parameters applied to this molecule
    for mol_idx, mol_forces in enumerate(molecule_force_list):
        # print(f'Forces for molecule {mol_idx}')
        for force_tag, force_dict in mol_forces.items():
            if force_tag == "ProperTorsions":
                for (atom_indices, parameter) in force_dict.items():
                    params.append(parameter.id)
                    indices.append(atom_indices)
                    #torsion_indices=tuple(torsion_indices)
                    #print(type(torsion_indices))
                    print(torsion_indices)
                    #print(type(atom_indices))
                    print(atom_indices)
                    if atom_indices == torsion_indices or tuple(
                        reversed(atom_indices)
                    ) == torsion_indices:
                        #mol.SetData("IDMatch", parameter.id)
                        tid=parameter.id
    print(params)
    print(indices)
    return tid


client = ptl.FractalClient()
# from the TorsionDriveDataset collection picking up given datasetName
ds = client.get_collection("TorsionDriveDataset", 'OpenFF Substituted Phenyl Set 1')

def testQuery(smiles):
    #print(ds.get_entry(smiles))
    #print(dir(ds.get_entry(smiles)))
    dih=ds.get_entry(smiles).dict()['td_keywords']['dihedrals'][0]
    print(dih)
    mapped_smiles = ds.get_entry(smiles).attributes["canonical_isomeric_explicit_hydrogen_mapped_smiles"]
    #print(mapped_smiles)
    return mapped_smiles, dih


def biphenyl(filename):
    with open(filename) as json_file:
        data = json.load(json_file)
    for key, item in data.items():
        testQuery(key)
biphenyl('biphenyls_set_input.json')


color_keys= ['maroon', 'brown', 'indianred', 'red', 'coral','orange', 'gold', 'darkkhaki','yellowgreen','limegreen',
              'mediumseagreen', 'teal', 'steelblue', 'cornflowerblue', 'royalblue', 'darkblue',
              'mediumblue', 'slateblue', 'blueviolet', 'purple','mediumvioletred', 'deeppink', 'hotpink',
              'palevioletred', 'pink', 'lightpink']



color_keys2=['darkblue',
              'mediumblue', 'slateblue', 'blueviolet', 'purple','mediumvioletred', 'deeppink', 'hotpink',
              'cornflowerblue', 'pink', 'lightpink']

color_keys2=['teal', 'hotpink', 'purple', 'gold', 'orange', 'slateblue', 'darkkhaki', 'lightpink', 'purple', 'hotpink']

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
    plt.ylabel(fgroup_symbols_colors[fgroup][0], rotation=0, size=10,
               color=colors[fgroup_symbols_colors[fgroup][1]])
    if i == len(fgroup_wbos)-1:
        plt.xlabel('AM1 ELF10 Wiberg bond order', fontsize=14)
        plt.xticks(fontsize=14)
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
fontsize=8
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
                     textcoords='offset points', xytext=(3, 2), color='black', fontsize=fontsize)
            continue
        plt.plot(x, y, '.', markersize=markersize, color=fgroup_symbols_colors[row[1].substituent][1],
                     label=fgroup_symbols_colors[row[1].substituent][0])
        plt.annotate(fgroup_symbols_colors[row[1].substituent][0], (x, y),
                     textcoords='offset points', xytext=(3, 2), color= fgroup_symbols_colors[row[1].substituent][1], fontsize=fontsize)

    plt.xlim(0.83, 1.16)
    plt.ylim(-0.86, 0.85)
    plt.ylabel(r'$\sigma_{}$'.format(sigma), fontsize=14)
    plt.xlabel('AM1 ELF10 Wiberg Bond Order', fontsize=14);
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    if sigma == 'm':
        r_value = df.corr().sigma_m.wbo_r_meta
    if sigma == 'p':
        r_value = df.corr().sigma_p.wbo_r_para
    #print(r_value)
    textstr = r'$\rho =%.2f$' % (r_value)
    props = dict(boxstyle='square', facecolor='white', alpha=0.5)
    ax.text(0.75, 0.95, textstr, transform=ax.transAxes, fontsize=14,
            verticalalignment='top', bbox=props)
    plt.tight_layout()
    fig.savefig('figures/hammett_sigma_{}.pdf'.format(sigma))


# Generate torsion barrier height vs ELF10 AM1 WBO plot
with open('../../phenyl_benchmark/data/qcarchive_torsiondrives.json', 'r') as f:
    fgroups_td = json.load(f)

# Generate 2 plots. One for good lines and one for lines that have issues
plot_1 = ['dimethylamino', 'methylamino', 'ethylamino', 'propylamino', 'hydroxy', 'methoxy', 'phenylurea', 'benzoicacid', 'nitro']
plot_2 = ['amino', 'ethoxy', 'dimethylurea', 'urea', 'ethylamide', 'amide', 'carbamate', 'ethoxycarbonyl']
symbols = ['o', 'P', '^', '*', 's', 'p',  'X', 'd', 'H', '>']

both_plots=plot_1 + plot_2

def r_value_ci(am1_wbos, max_energies):
    return stats.linregress(am1_wbos, max_energies)[2]**2

fontsize = 14
fig, ax = plt.subplots()
colors = []
r_values = []
for i, fgroup in enumerate(plot_1):
    if fgroup not in fgroups_td:
        print(fgroup)
        continue
    energies = fgroups_td[fgroup]['energy']
    am1_wbos = fgroups_td[fgroup]['elf10_am1_wbo']
    max_energies = [max(energy) for energy in energies]
    slope, intercept, r_value, p_value, std_err = stats.linregress(am1_wbos, max_energies)
    r_ci = arch.bootstrap.IIDBootstrap(np.asarray(am1_wbos), np.asarray(max_energies)).conf_int(r_value_ci, 1000, method='percentile')
    #print(r_ci)
    fgroups_td[fgroup]['stats'] = [slope, std_err, r_value**2, r_ci[0][0], r_ci[1][0]]
    plt.plot(np.unique(am1_wbos), np.poly1d([slope, intercept])(np.unique(am1_wbos)), fgroup_symbols_colors[fgroup][1])
    plt.scatter(x=am1_wbos, y=max_energies, color=fgroup_symbols_colors[fgroup][1], marker=symbols[i], label=fgroup_symbols_colors[fgroup][0])
    colors.append(fgroup_symbols_colors[fgroup][1])
    r_values.append([r_value**2, r_ci[0][0], r_ci[1][0]])

l = ax.legend(bbox_to_anchor=(1, 1), fontsize=fontsize)
for i, text in enumerate(l.get_texts()):
    text.set_color(colors[i])

plt.xlabel('AM1 ELF10 Wiberg bond order', fontsize=fontsize)
plt.ylabel('Torsion barrier height (kJ/mol)', fontsize=fontsize)
plt.xlim(0.8, 1.3)
plt.ylim(0, 50)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.tight_layout()
plt.savefig('figures/energy_vs_wbo_1.pdf')

colors = []
ig, ax = plt.subplots()
tig_dict={'TIG1':[[],[]], 'TIG2':[[],[]], 'TIG3':[[],[]], 'TIG4':[[],[]], 'TIG5':[[],[]], 'TIG6':[[],[]], 'TIG7':[[],[]], 'TIG8':[[],[]], 'TIG9':[[],[]], 'TIG10':[[],[]]}
molDict={}
"""
for i, fgroup in enumerate(both_plots):
    if fgroup not in fgroups_td:
        continue
    print(i)
    print(fgroup)
    energies = fgroups_td[fgroup]['energy']
    am1_wbos = fgroups_td[fgroup]['elf10_am1_wbo']
    max_energies = [max(energy) for energy in energies]
    molcount=0
    torsions=[]
    for i, smiles in enumerate(fgroups_td[fgroup]['indices']):
        molDict[smiles]=[am1_wbos[i], max_energies[i]]
        molcount+=1
        #testQuery(smiles)
        #with open('../../phenyl_benchmark/data/{}_td_job_indices.json'.format(fgroup), 'r') as f:
        #/Users/jessica/Documents/Grad_research/fragmenter_data/wbo-manuscript-figures/proof_of_concept/data/data
        with open('data/data/{}_td_job_indices.json'.format(fgroup), 'r') as f:
            indices = json.load(f)
            for m in indices:
                if m[0] == smiles:
                    molDict[smiles].extend([m[1], m[4]])
    for sm, dd in molDict.items():
        print(dd)
        smiles, dih=testQuery(sm)
        ff='tig_proof_of_concept_1.3.0.offxml'
        tid = checkTorsion(smiles, dih, ff)
        torsions.append(tid)
        tig_dict[tid][0].append(dd[0])
        tig_dict[tid][1].append(dd[1])
    print(molcount)
    print(tig_dict)
    print(torsions)
    print(len(torsions))
    with open('biphenyl_data.pickle', 'rb') as handle:
        b = pickle.load(handle)
    for key, item in b.items():
        smiles, dih=testQuery(key)
        tid = checkTorsion(smiles, item[2], ff)
        tig_dict[tid][0].append(item[0])
        tig_dict[tid][1].append(item[1])

    import pickle
    with open("wbotb.pkl", "wb") as f:
        pickle.dump(tig_dict, f)
"""
def makeCovPlot(filename):
    with open(filename, "rb") as f:
        plotdata = pickle.load(f)
    #print(plotdata)
    count=0
    colors=[]
    tid_td={}
    for key, data in plotdata.items():
        am1_wbos=data[0]
        max_energies=data[1]
        if am1_wbos==[]:
            continue
        #print(am1_wbos)
        #print(max_energies)

        slope, intercept, r_value, p_value, std_err = stats.linregress(am1_wbos, max_energies)
        r_ci = arch.bootstrap.IIDBootstrap(np.asarray(am1_wbos), np.asarray(max_energies)).conf_int(r_value_ci, 10000, method='percentile')
        #print(r_ci)
        fgroups_td[fgroup]['stats'] = [slope, std_err, r_value**2, r_ci[0][0], r_ci[1][0]]
        tid_td[key] = [slope, std_err, r_value**2, r_ci[0][0], r_ci[1][0]]
        plt.plot(np.unique(am1_wbos), np.poly1d([slope, intercept])(np.unique(am1_wbos)), color_keys2[count])
        plt.scatter(x=am1_wbos, y=max_energies, color=color_keys2[count], marker=symbols[count], label=key)
        colors.append(color_keys2[count])
        count+=1
    #store statistics from the td vs wbo plot for table generation
    with open("table_data.pkl", "wb") as f:
        pickle.dump(tid_td, f)

    l = ax.legend(bbox_to_anchor=(1, 1), fontsize=fontsize)
    for i, text in enumerate(l.get_texts()):
        text.set_color(colors[i])

    plt.xlabel('AM1 ELF10 Wiberg bond order', fontsize=fontsize)
    plt.ylabel('Torsion barrier height (kJ/mol)', fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    #plt.xlim(0.8, 1.5)
    #plt.ylim(0, 100)
    plt.tight_layout()
    plt.savefig('energy_vs_wbo_full_newcolors.pdf')

makeCovPlot('wbotb.pkl')


# generate table
stats_table = {'Parameter': [], 'smarts':[], 'slope': [],'standard error': [],  'r^2': [], 'CI_1': [], 'CI_2': []}
#[slope, std_err, r_value**2, r_ci[0][0], r_ci[1][0]]

with open('table_data.pkl', 'rb') as f:
    tabledata = pickle.load(f)
smartsDict={
        'TIG1':'[*:1]~[#6X3:2]-[#6X3:3]~[*:4]',
        'TIG2':'[*:1]~[#6X3:2]-[#6X3$(*=[#8,#16,#7]):3]~[*:4]',
        'TIG3':'[*:1]~[#6X3:2]-[#6X3:3](-[#8H1])=[#8X1:4]',
        'TIG4':'[*:1]~[#7X3:2]-!@[#6X3:3]~@[#6:4]',
        'TIG5':'[#6X3:1]~[#7X3:2]-!@[#6X3:3]~@[#6:4]',
        'TIG6':'[#6X3$(*~[#6]):1]~[#7X3:2]-!@[#6X3:3]~@[#6:4]',
        'TIG7':'[#6X4:1]~[#7X3:2]-!@[#6X3:3]~@[#6:4]',
        'TIG8':'[#8X1:1]~[#7X3:2]~[#6X3:3]~[*:4]',
        'TIG9':'[*:1]~[#6X3:2]-[#8X2:3]-[*:4]',
        'TIG10':'[*:1]~[#6X3:2]-[#8X2:3]-[#1:4]'
        }
for key, item in tabledata.items():
    stats_table['Parameter'].append(key)
    stats_table['smarts'].append(smartsDict[key])
    stats_table['slope'].append(round(item[0],2))
    stats_table['standard error'].append(round(item[1],2))
    stats_table['r^2'].append(round(item[2],2))
    stats_table['CI_1'].append(round(item[3], 2))
    stats_table['CI_2'].append(round(item[4], 2))
latex_table = pd.DataFrame(stats_table).to_latex(index=False)
with open('figures/stats_tid.tex', 'w') as f:
    f.write(latex_table)

