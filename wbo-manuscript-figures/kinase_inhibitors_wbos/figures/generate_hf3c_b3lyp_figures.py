import qcfractal.interface as ptl
from openeye import oechem
import fragmenter
import cmiles

import numpy as np
import json
import argparse

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sbn
import arch.bootstrap

# Plotting functions
def get_atom_map(molecule, mapped_smiles):
    """
    The canonical mapping on oemols and cmiles mapped SMILES does not lend itself to easy sorting of bonds
    for figures. Some old tagged SMILES I have (generated from SMILES that were not canonical but in an order
    that reads easily left to right) have the atoms mapped in such a way that it is easy to sort the bonds for
    figures. This function is taken from cmiles with a small change. The map goes from {map_in_mol: map_in_tagged_smiles}

    Parameters
    ----------
    molecule: oemol with map indices (tags were added using canonical ordering from cmiles)
    mapped_smiles: tagged SMILES to map atom maps to

    Returns
    -------
    atom_map: dict
        {mol_map_idx: smiles_map}
    """
    # check that smiles has explicit hydrogen and map indices
    mapped_mol = oechem.OEMol()
    oechem.OESmilesToMol(mapped_mol, mapped_smiles)
    if not cmiles.utils.has_atom_map(mapped_mol):
        raise ValueError("Mapped SMILES must have map indices for all atoms and hydrogens")
    # Check molecule for explicit hydrogen
    if not cmiles.utils.has_explicit_hydrogen(molecule) and strict:
        raise ValueError("Molecule must have explicit hydrogens")

    # canonical order mapped mol to ensure atom map is always generated in the same order
    cmiles._cmiles_oe.canonical_order_atoms(mapped_mol)
    aopts = oechem.OEExprOpts_DefaultAtoms
    bopts = oechem.OEExprOpts_DefaultBonds
    ss = oechem.OESubSearch(mapped_mol, aopts, bopts)
    oechem.OEPrepareSearch(molecule, ss)
    ss.SetMaxMatches(1)

    atom_map = {}
    matches = [m for m in ss.Match(molecule)]
    if not matches:
        raise RuntimeError("MCSS failed for {}, smiles: {}".format(oechem.OEMolToSmiles(molecule), mapped_smiles))
    for match in matches:
        for ma in match.GetAtoms():
            atom_map[ma.target.GetMapIdx()] = ma.pattern.GetMapIdx()

    # sanity check
    mol = oechem.OEGraphMol()
    oechem.OESubsetMol(mol, match, True)
    matched_smiles = cmiles.utils.mol_to_smiles(mol, isomeric=False, explicit_hydrogen=False, mapped=False)
    molcopy = oechem.OEMol(molecule)
    smiles = cmiles.utils.mol_to_smiles(molcopy, isomeric=False, explicit_hydrogen=False, mapped=False)
    pattern_smiles = cmiles.utils.mol_to_smiles(mapped_mol, isomeric=False, explicit_hydrogen=False, mapped=False)
    if not matched_smiles == smiles == pattern_smiles:
        raise RuntimeError("Matched molecule, input molecule and mapped SMILES are not the same ")
    return atom_map

def plot_distributions(bond_keys, bond_wbos, bond_order_type, spec, fname):
    """
    Plot distribution of wbos
    bond_keys: list of tuples for bonds to plot
    bond_wbos: dict of bond_tuple:wbos (return from get_bond_wbos)
    bond_order_type: str
        wiberg_lowdin or mayer
    name: str, name of molecule
    """
    # Generate WBO distributions
    colors = fragmenter.chemi._KELLYS_COLORS
    n = len(bond_keys)
    fig, axes = plt.subplots(n, 1)
    fig.dpi = 400
    x_min = 3
    x_max = 0
    sorted_bond_keys = sorted(list(bond_keys.keys()))
    for b in sorted_bond_keys:
        try:
            wbo = bond_wbos['rotors'][bond_keys[b]][bond_order_type]
        except KeyError:
            wbo = bond_wbos['others'][bond_keys[b]][bond_order_type]
        if min(wbo) < x_min:
            x_min = min(wbo)
        if max(wbo) > x_max:
            x_max = max(wbo)

    # Format labels
    f = mticker.ScalarFormatter(useOffset=False, useMathText=True)
    g = lambda x,pos : "${}$".format(f._formatSciNotation('%1.10e' % x))
    fmt = mticker.FuncFormatter(g)
    for i, bond in enumerate(sorted_bond_keys):
        try:
            wbo = bond_wbos['rotors'][bond_keys[bond]][bond_order_type]
        except KeyError:
            wbo = bond_wbos['others'][bond_keys[bond]][bond_order_type]
        var = np.var(wbo)
        ci = arch.bootstrap.IIDBootstrap(np.asarray(wbo)).conf_int(np.var, 1000)
        ax = plt.subplot(n, 1, i+1)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.patch.set_facecolor('none')
        # format values
        var_txt = '{0:.2g}'.format(var/1e-5)
        ci_0_txt = '{0:.2g}'.format(ci[0][0]/1e-5)
        ci_1_txt = '{0:.2g}'.format(ci[1][0]/1e-5)
        textstr = r'%s$_{%s}^{%s}$' % (var_txt, ci_0_txt, ci_1_txt)
        if i==0:
            ax.text(-0.11, 0.84, 'Bonds   Variance (1E-5)', transform=ax.transAxes, fontsize=12,
                    verticalalignment='top')#, bbox=props)
        sbn.kdeplot(wbo, shade= True, alpha=0.85, color=colors[i], label=textstr)
        l = plt.legend(handlelength=0, fontsize=12, loc='lower left', frameon=False, markerscale=0.1)
        plt.setp(l.get_texts(), color=colors[i])
        plt.xlim(0.65, 1.45)
        plt.xticks(fontsize=14)
        plt.yticks([])
        ax.yaxis.set_label_coords(-0.05, 0)
        plt.ylabel(bond, rotation=0, size=14, color=colors[i])
        if i != n-1:
            plt.xticks([])
        else:
            if bond_order_type == 'wiberg_lowdin':
                bond_order_type_str = 'Wiberg'
            else:
                bond_order_type_str = 'Mayer'
            plt.xlabel('{} bond order'.format(bond_order_type_str), fontsize=14)
        if i == 0:
            if bond_order_type == 'wiberg_lowdin':
                bond_order_type_str = 'Wiberg'
            else:
                bond_order_type_str = 'Mayer'
            #plt.legend(prop={'size': 10}, bbox_to_anchor=(1.35, 1))
            #plt.title("Wiberg Löwding distributions over conformations for highlighted bonds", fontsize=14)
            if spec == 'default':
                spec = 'b3lyp'
            plt.title('{} {} bond order distributions'.format(spec.upper(), bond_order_type_str))

            overlap=0.5
    h_pad = 5 + (- 5*(1 + overlap))
    fig.tight_layout(h_pad=h_pad)
    plt.savefig(fname)

def sort_rings(tagged_smiles):
    """
    Sort ring bonds for nicer plottin
    Paramerts
    ---------
    tagged_smiles: mapped SMILES for plotting

    Returns
    -------
    ring_bonds: list of tuples of aliphatic and aromtic ring bonds atom map indices
    ring_sizes: list of ring sizes
    """
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, tagged_smiles)
    nringsystems, parts = oechem.OEDetermineRingSystems(mol)
    aromatic_ring_bonds = []
    aliphatic_ring_bonds = []
    aromatic_ring_sizes = []
    aliphatic_ring_sizes = []
    for ring_idx in range(1, nringsystems + 1):
        i = 0
        j = 0
        for bond in mol.GetBonds():
            if bond.IsInRing():
                if bond.IsAromatic():
                    if parts[bond.GetBgnIdx()] == ring_idx and parts[bond.GetEndIdx()] == ring_idx:
                        aromatic_ring_bonds.append((bond.GetBgn().GetMapIdx(), bond.GetEnd().GetMapIdx()))
                        i += 1

                if not bond.IsAromatic():
                    if parts[bond.GetBgnIdx()] == ring_idx and parts[bond.GetEndIdx()] == ring_idx:
                        aliphatic_ring_bonds.append((bond.GetBgn().GetMapIdx(), bond.GetEnd().GetMapIdx()))
                        j += 1
        if i > 0:
            aromatic_ring_sizes.append(i)
        if j > 0:
            aliphatic_ring_sizes.append(j)

    return aliphatic_ring_bonds + aromatic_ring_bonds, aliphatic_ring_sizes + aromatic_ring_sizes

def correlation_plot(sorted_bonds, bonds_mappings, bonds_wbos, bo_type, spec, ring_sizes, fname):
    """
    Generate correlation plot of bond WBOs
    parameters
    ----------
    sorted_bonds: list of bond tuples sorted for plotting
    bonds_mapping: dict
        map of bond tuple in sorted bonds to corresponding bond tuple in molecule
    bonds_wbos: dict
        mapping of bond tuple (in molecule) to list of WBOs
    bo_type: str
     wiberg_lowdin or mayer
    fname: str
        filename
    """
    correlation_coef = np.zeros((len(sorted_bonds), len(sorted_bonds)))

    for i, bond_1 in enumerate(sorted_bonds):
        for j, bond_2 in enumerate(sorted_bonds):
            try:
                bond_1_mapping = bonds_mappings[bond_1]
            except KeyError:
                bond_1_mapping = bonds_mappings[tuple(reversed(bond_1))]
            try:
                bond_2_mapping = bonds_mappings[bond_2]
            except KeyError:
                bond_2_mapping = bonds_mappings[tuple(reversed(bond_2))]
            bond_1_wbo = bonds_wbos[bond_1_mapping][bo_type]
            bond_2_wbo = bonds_wbos[bond_2_mapping][bo_type]
            corr_coef = np.corrcoef(bond_1_wbo, bond_2_wbo)
            correlation_coef[i][j] = corr_coef[0, 1]

    fig, ax = plt.subplots()
    corr = ax.imshow(correlation_coef, cmap='coolwarm')
    ax.xaxis.set_ticks_position('bottom')

    plt.xticks(np.arange(len(sorted_bonds)), sorted_bonds, rotation='vertical', size=8);
    plt.yticks(np.arange(len(sorted_bonds)), sorted_bonds, size=8);

    linewidth = 1.0
    x_1 = len(sorted_bonds) - np.cumsum(ring_sizes)[-1] - 0.5

    ax.axvline(x=x_1, color='white', linewidth=linewidth)
    ax.axhline(y=x_1, color='white', linewidth=linewidth)
    for i in ring_sizes:
        x_1 += i
        ax.axvline(x=x_1, color='white', linewidth=linewidth)
        ax.axhline(y=x_1, color='white', linewidth=linewidth)

    fig.colorbar(corr)
    if spec == 'default':
        spec = 'B3LYP'
    else:
        spec = 'HF3C'
    if bo_type == 'wiberg_lowdin':
        bo_type = 'WBO'
    else:
        bo_type = 'MBO'
    ax.set_title('{} {} correlations'.format(spec, bo_type));
    fig.savefig(fname)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run psi4 to calculate bond orders')
    parser.add_argument('-n', '--name', type=str,
                        help='kinase inhibitor')

    args = parser.parse_args()
    name = args.name

    # Get old tagged SMILES. These were generated from non canonical SMILES (they were taken from drug bank) where
    # the molecules created from these SMILES had the atoms in a nice order for visualization and plotting. Canonical
    # order (from cmiles) is not that pleasing for visualization so here I use the tagged SMILES to map the data onto
    with open('../data_generation/archive/clinical_kinase_inhibitors_tagged.smi', 'r') as f:
        smiles = f.read()
    smiles = smiles.split('\n')
    for sm in smiles:
        sm = sm.split(' ')
        if sm[-1] == name:
            tagged_smiles = sm[0]
    # Get dataset
    client = ptl.FractalClient()
    ds = client.get_collection('OptimizationDataset', 'Kinase Inhibitors: WBO Distributions')

    with open('../data_generation/data/kinase_inhibitors_cmiles_ids.json', 'r') as f:
        identifiers = json.load(f)
    cmiles_ids = identifiers[name]
    idx_smiles = cmiles_ids['canonical_isomeric_smiles']

    # # Download all Bond orders and save
    # bond_wbos = {'hf3c': {'rotors': {}, 'rings': {}, 'others': {}},
    #              'default': {'rotors': {}, 'rings': {}, 'others': {}}}
    # for index in ds.df.index:
    #     if idx_smiles in index:
    #         print(index)
    #         entry = ds.get_entry(index)
    #         mapped_smiles = entry.attributes['canonical_isomeric_explicit_hydrogen_mapped_smiles']
    #         mapped_mol = oechem.OEMol()
    #         oechem.OESmilesToMol(mapped_mol, mapped_smiles)
    #         n_atoms = mapped_mol.GetMaxAtomIdx()
    #         for spec in bond_wbos:
    #             rec = ds.get_record(index, spec)
    #             opt = rec.get_trajectory()[-1]
    #             wiberg = np.array(opt.extras['qcvars']['WIBERG_LOWDIN_INDICES']).reshape(-1, n_atoms)
    #             mayer = np.array(opt.extras['qcvars']['MAYER_INDICES']).reshape(-1, n_atoms)
    #             for bond in mapped_mol.GetBonds():
    #                 a_1 = bond.GetBgn()
    #                 a_2 = bond.GetEnd()
    #                 if a_1.IsHydrogen() or a_2.IsHydrogen():
    #                     continue
    #                 idx_1 = a_1.GetMapIdx()
    #                 idx_2 = a_2.GetMapIdx()
    #                 wbo = wiberg[idx_1 - 1][idx_2 - 1]
    #                 mbo = mayer[idx_1 - 1][idx_2 - 1]
    #                 bond_tuple = (idx_1, idx_2)
    #                 if bond.IsRotor():
    #                     if bond_tuple not in bond_wbos[spec]['rotors']:
    #                         bond_wbos[spec]['rotors'][bond_tuple] = {'wiberg_lowdin': [], 'mayer': []}
    #                     bond_wbos[spec]['rotors'][bond_tuple]['wiberg_lowdin'].append(wbo)
    #                     bond_wbos[spec]['rotors'][bond_tuple]['mayer'].append(mbo)
    #                 elif bond.IsInRing():
    #                     if bond_tuple not in bond_wbos[spec]['rings']:
    #                         bond_wbos[spec]['rings'][bond_tuple] = {'wiberg_lowdin': [], 'mayer': []}
    #                     bond_wbos[spec]['rings'][bond_tuple]['wiberg_lowdin'].append(wbo)
    #                     bond_wbos[spec]['rings'][bond_tuple]['mayer'].append(mbo)
    #                 else:
    #                     if bond_tuple not in bond_wbos[spec]['others']:
    #                         bond_wbos[spec]['others'][bond_tuple] = {'wiberg_lowdin': [], 'mayer': []}
    #                     bond_wbos[spec]['others'][bond_tuple]['wiberg_lowdin'].append(wbo)
    #                     bond_wbos[spec]['others'][bond_tuple]['mayer'].append(mbo)
    # # serialize and save
    # bond_wbos_ser = {}
    # for spec in bond_wbos:
    #     bond_wbos_ser[spec] = {}
    #     for bond_type in bond_wbos[spec]:
    #         bond_wbos_ser[spec][bond_type] = {}
    #         for bond_tuple in bond_wbos[spec][bond_type]:
    #             ser_bond = fragmenter.utils.serialize_bond(bond_tuple)
    #             bond_wbos_ser[spec][bond_type][ser_bond] = bond_wbos[spec][bond_type][bond_tuple]
    # fname = '../data_generation/data/{}_hf3c_b3lyp_wbos.json'.format(name)
    # with open(fname, 'w') as f:
    #     json.dump(bond_wbos_ser, f, indent=2, sort_keys=True)

    mapped_smiles = cmiles_ids['canonical_isomeric_explicit_hydrogen_mapped_smiles']
    mapped_mol = oechem.OEMol()
    oechem.OESmilesToMol(mapped_mol, mapped_smiles)

    with open('../data_generation/data/{}_hf3c_b3lyp_wbos.json'.format(name), 'r') as f:
        bond_wbos_ser = json.load(f)
    bond_wbos = {}
    for spec in bond_wbos_ser:
        bond_wbos[spec] = {}
        for bond_type in bond_wbos_ser[spec]:
            bond_wbos[spec][bond_type] = {}
            for bond_ser in bond_wbos_ser[spec][bond_type]:
                bond_tuple = fragmenter.utils.deserialize_bond(bond_ser)
                bond_wbos[spec][bond_type][bond_tuple] = bond_wbos_ser[spec][bond_type][bond_ser]


    # Some data munging for plotting
    # First get old atom map that is more pleasing for visualization
    with open('../data_generation/archive/clinical_kinase_inhibitors_tagged.smi', 'r') as f:
        smiles = f.read()
    smiles = smiles.split('\n')
    for sm in smiles:
        sm = sm.split(' ')
        if sm[-1] == name:
            old_mapped_smiles = sm[0]

    atom_map = get_atom_map(mapped_mol, old_mapped_smiles)
    for spec in bond_wbos:
        if spec == 'default':
            dir_name = 'b3lyp'
        else:
            dir_name = spec
        mapped_rotor_bonds = {(atom_map[bond[0]], atom_map[bond[1]]): bond for bond in bond_wbos[spec]['rotors']}
        mapped_other_bonds = {(atom_map[bond[0]], atom_map[bond[1]]): bond for bond in bond_wbos[spec]['others']}
        mapped_ring_bonds = {(atom_map[bond[0]], atom_map[bond[1]]): bond for bond in bond_wbos[spec]['rings']}

        bonds_to_highlight = {**mapped_rotor_bonds, **mapped_other_bonds}
        bonds_to_highlight_sorted = sorted(list(bonds_to_highlight.keys()))
        rotor_bonds_sorted = sorted(list(mapped_rotor_bonds.keys()))

        fragmenter.chemi.highlight_bond_by_map_idx(old_mapped_smiles, bond_map_idx=bonds_to_highlight_sorted,
                                                   map_idx=True,
                                                   label_scale=1.5,
                                                   fname='{}_figures/{}_all_single_bonds_highlighted.pdf'.format(dir_name, name))
        fragmenter.chemi.highlight_bond_by_map_idx(old_mapped_smiles, bond_map_idx=rotor_bonds_sorted, map_idx=True,
                                                   label_scale=1.5, fname='{}_figures/{}_rotor_bonds_highlighted.pdf'.format(dir_name, name))

        for bo_type in ('wiberg_lowdin', 'mayer'):
            plot_distributions(bonds_to_highlight, bond_wbos[spec], spec=spec,
                               bond_order_type=bo_type, fname='{}_figures/{}_all_single_{}_distributions.pdf'.format(dir_name, name, bo_type))
            plot_distributions(mapped_rotor_bonds, bond_wbos[spec], spec=spec,
                               bond_order_type=bo_type, fname='{}_figures/{}_rotor_{}_distributions.pdf'.format(dir_name, name, bo_type))

        # More data munging to get bonds sorted in order for correlation plot
        sorted_ring_bonds, ring_sizes = sort_rings(old_mapped_smiles)
        bonds_wbos_all = {**bond_wbos[spec]['rotors'], **bond_wbos[spec]['others'], **bond_wbos[spec]['rings']}
        sorted_bonds = bonds_to_highlight_sorted + sorted_ring_bonds
        bonds_mappings = {**bonds_to_highlight, **mapped_ring_bonds}

        for bo_type in ('wiberg_lowdin', 'mayer'):
            correlation_plot(sorted_bonds, bonds_mappings, bonds_wbos_all, bo_type, spec, ring_sizes,
                             fname='{}_figures/{}_{}_correlations.pdf'.format(dir_name, name, bo_type))
