import json
from openeye import oechem, oequacpac
import fragmenter
import cmiles
import glob
from openmoltools import openeye
import matplotlib.pyplot as plt
import seaborn as sbn
import numpy as np
import argparse

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

def get_bond_wbos(mols_wbos):
    """
    Collect all bonds wbos from list of oemols that have wbos calculated

    Parameters
    ----------
    mols_wbos: list of mapped oemols

    Returns
    -------
    bonds: dict
        maps type of bonds --> bonds
                bonds --> list of WBOs
        bond --> (map_idx_1, map_idx_2)
    """
    bonds = {'rotors': {}, 'rings': {}, 'others': {}}
    for mol in mols_wbos:
        for bond in mol.GetBonds():
            a1 = bond.GetBgn()
            a2 = bond.GetEnd()
            if not a1.IsHydrogen() and not a2.IsHydrogen():
                if bond.IsRotor():
                    key = (bond.GetBgn().GetMapIdx(), bond.GetEnd().GetMapIdx())
                    if key not in bonds['rotors']:
                        bonds['rotors'][key] = []
                    bonds['rotors'][key].append(bond.GetData('WibergBondOrder'))
                elif bond.IsInRing():
                    key = (bond.GetBgn().GetMapIdx(), bond.GetEnd().GetMapIdx())
                    if key not in bonds['rings']:
                        bonds['rings'][key] = []
                    bonds['rings'][key].append(bond.GetData('WibergBondOrder'))
                else:
                    key = (bond.GetBgn().GetMapIdx(), bond.GetEnd().GetMapIdx())
                    if key not in bonds['others']:
                        bonds['others'][key] = []
                    bonds['others'][key].append(bond.GetData('WibergBondOrder'))
    return bonds

def plot_distributions(bond_keys, bond_wbos, name, fname):
    """
    Plot distribution of wbos
    bond_keys: list of tuples for bonds to plot
    bond_wbos: dict of bond_tuple:wbos (return from get_bond_wbos)
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
            wbo = bond_wbos['rotors'][bond_keys[b]]
        except KeyError:
            wbo = bond_wbos['others'][bond_keys[b]]
        if min(wbo) < x_min:
            x_min = min(wbo)
        if max(wbo) > x_max:
            x_max = max(wbo)

    #sorted_rot_bonds = sorted(bonds['rotors'].keys())
    for i, bond in enumerate(sorted_bond_keys):
        try:
            wbo = bond_wbos['rotors'][bond_keys[bond]]
        except KeyError:
            wbo = bond_wbos['others'][bond_keys[bond]]
        ax = plt.subplot(n, 1, i+1)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.patch.set_facecolor('none')
        sbn.kdeplot(wbo, shade= True, alpha=0.85, color=colors[i])
      #sbn.distplot(wbo, hist=False, rug=True, kde=False, color='black')
        sbn.kdeplot(wbo, lw=1.5, color=colors[i])
        #plt.axvline(x=wbo_s, ymin=0, ymax=1, color='black', linewidth=0.5)
        plt.xlim(x_min-0.05, x_max+0.05)
        plt.yticks([])
        #ax.set_yticklabels(bond_order_std_rot_bonds_opt['Imatinib']['bonds'], fontsize=6, rotation=0)
        #ax.yaxis.grid(False)
        ax.yaxis.set_label_coords(-0.05, 0)
        plt.ylabel(bond, rotation=0, size=8)
        if i != n-1:
            plt.xticks([])
        else:
            plt.xlabel('Bond order')
        if i == 0:
            #plt.legend(prop={'size': 10}, bbox_to_anchor=(1.35, 1))
            plt.title("{} WBO".format(name))
            overlap=0.5
    h_pad = 5 + (- 5*(1 + overlap))
    fig.tight_layout(h_pad=h_pad)
    plt.savefig(fname)

def sort_rings(tagged_smiles):
    """
    Sort ring bonds for nicer plotting

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

def correlation_plot(sorted_bonds, bonds_mapping, bonds_wbos, ring_sizes, fname):
    """
    Generate correlation plot of bond WBOs
    parameters
    ----------
    sorted_bonds: list of bond tuples sorted for plotting
    bonds_mapping: dict
        map of bond tuple in sorted bonds to corresponding bond tuple in molecule
    bonds_wbos: dict
        mapping of bond tuple (in molecule) to list of WBOs
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
            bond_1_wbo = bonds_wbos[bond_1_mapping]
            bond_2_wbo = bonds_wbos[bond_2_mapping]
            corr_coef = np.corrcoef(bond_1_wbo, bond_2_wbo)
            correlation_coef[i][j] = corr_coef[0, 1]

    fig, ax = plt.subplots()
    corr = ax.imshow(correlation_coef, cmap='coolwarm')
    ax.xaxis.set_ticks_position('bottom')

    plt.xticks(np.arange(len(sorted_bonds)), sorted_bonds, rotation='vertical', size=5.0);
    plt.yticks(np.arange(len(sorted_bonds)), sorted_bonds, size=5.0);

    linewidth = 1.0
    x_1 = len(sorted_bonds) - np.cumsum(ring_sizes)[-1] - 0.5

    ax.axvline(x=x_1, color='white', linewidth=linewidth)
    ax.axhline(y=x_1, color='white', linewidth=linewidth)
    for i in ring_sizes:
        x_1 += i
        ax.axvline(x=x_1, color='white', linewidth=linewidth)
        ax.axhline(y=x_1, color='white', linewidth=linewidth)

    fig.colorbar(corr)
    ax.set_title('{} {} conformers'.format(name, len(bond_1_wbo)));
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

    oemols = fragmenter.chemi.file_to_oemols('../data_generation/data/{}_am1_wbo.oeb'.format(name))
    atom_map = get_atom_map(oemols[0], tagged_smiles)

    bonds_wbos = get_bond_wbos(oemols)

    # Some data munging to get bonds sorted in the right order for better plotting
    mapped_rotor_bonds = {(atom_map[bond[0]], atom_map[bond[1]]): bond for bond in bonds_wbos['rotors']}
    mapped_other_bonds = {(atom_map[bond[0]], atom_map[bond[1]]): bond for bond in bonds_wbos['others']}
    mapped_ring_bonds = {(atom_map[bond[0]], atom_map[bond[1]]): bond for bond in bonds_wbos['rings']}

    bonds_to_highlight = {**mapped_rotor_bonds, **mapped_other_bonds}
    bonds_to_highlight_sorted = sorted(list(bonds_to_highlight.keys()))
    rotor_bonds_sorted = sorted(list(mapped_rotor_bonds.keys()))

    fragmenter.chemi.highlight_bond_by_map_idx(tagged_smiles, bond_map_idx=bonds_to_highlight_sorted, map_idx=True,
                                           label_scale=1.5, fname='am1_figures/{}_all_single_bonds_highlighted.pdf'.format(name))
    fragmenter.chemi.highlight_bond_by_map_idx(tagged_smiles, bond_map_idx=rotor_bonds_sorted, map_idx=True,
                                           label_scale=1.5, fname='am1_figures/{}_rotor_bonds_highlighted.pdf'.format(name))

    plot_distributions(bonds_to_highlight, bonds_wbos, name, fname='am1_figures/{}_all_single_wbo_distributions.pdf'.format(name))
    plot_distributions(mapped_rotor_bonds, bonds_wbos, name, fname='am1_figures/{}_rotor_wbo_distributions.pdf'.format(name))


    # More data munging to get bonds sorted in order for correlation plot
    sorted_ring_bonds, ring_sizes = sort_rings(tagged_smiles)
    bonds_wbos_all = {**bonds_wbos['rotors'], **bonds_wbos['others'], **bonds_wbos['rings']}
    sorted_bonds = bonds_to_highlight_sorted + sorted_ring_bonds
    bonds_mappings = {**bonds_to_highlight, **mapped_ring_bonds}

    correlation_plot(sorted_bonds, bonds_mappings, bonds_wbos_all, ring_sizes, fname='am1_figures/{}_wbo_correlations.pdf'.format(name))