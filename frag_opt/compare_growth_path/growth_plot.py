from openeye import oechem
from fragmenter import chemi, fragment, workflow_api
import cmiles._cmiles_oe as utils
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import copy
import json
import os

def sort_by_wbo(molecule, atoms, bonds):
    sort_by = []
    atoms_to_add = []
    for m_idx in atoms:
        a = molecule.GetAtom(oechem.OEHasMapIdx(m_idx))
        for nbr in a.GetAtoms():
            nbr_m_idx = nbr.GetMapIdx()
            if not nbr.IsHydrogen() and not nbr_m_idx in atoms:
                b = molecule.GetBond(a, nbr)
                atoms_to_add.append((nbr_m_idx, (m_idx, nbr_m_idx)))
                sort_by.append(b.GetData('WibergBondOrder'))
    return atoms_to_add, sort_by


def add_one_bond(molecule, atoms, bonds, functional_groups, ring_systems, target_bond, heuristic):
    bond_atom_1 = molecule.GetAtom(oechem.OEHasMapIdx(target_bond[0]))
    bond_atom_2 = molecule.GetAtom(oechem.OEHasMapIdx(target_bond[1]))
    atoms_to_add = []
    sort_by_1 = []
    sort_by_2 = []
    sort_by = []
    for m_idx in atoms:
        a = molecule.GetAtom(oechem.OEHasMapIdx(m_idx))
        if a is None:
            print('{} is none'.format(m_idx))
            continue
        for nbr in a.GetAtoms():
            nbr_m_idx = nbr.GetMapIdx()
            if not nbr.IsHydrogen() and not nbr_m_idx in atoms:
                b = molecule.GetBond(a, nbr)
                atoms_to_add.append((nbr_m_idx, (m_idx, nbr_m_idx)))
                if heuristic == 'wbo':
                    sort_by.append(b.GetData('WibergBondOrder'))
                    reverse = True
                elif heuristic == 'path_length':
                    sort_by_1.append(oechem.OEGetPathLength(bond_atom_1, nbr))
                    sort_by_2.append(oechem.OEGetPathLength(bond_atom_2, nbr))
                    #sort_by.append(oechem.OEGetPathLength(bond_atom, nbr))
                    reverse = False
                else:
                    raise ValueError('Only wbo and path_lenght are supported heuristics')
    if heuristic == 'path_length':
        min_1 = min(sort_by_1)
        min_2 = min(sort_by_2)
        if min_1 < min_2:
            sort_by = sort_by_1
        elif min_2 < min_1:
            sort_by = sort_by_2
        elif min_1 == min_2:
            # Find the equivalent atoms and take the one with higher WBO
            indices = []
            for l in [sort_by_1, sort_by_2]:
                indices.extend([i for i, x in enumerate(l) if x == min_1])
            atoms_to_add = [atoms_to_add[i] for i in indices]
            for a, b in atoms_to_add:
                a1 = molecule.GetAtom(oechem.OEHasMapIdx(b[0]))
                a2 = molecule.GetAtom(oechem.OEHasMapIdx(b[1]))
                bond = molecule.GetBond(a1, a2)
                sort_by.append(bond.GetData('WibergBondOrder'))
                reverse = True

    sorted_atoms = [a for _, a in sorted(zip(sort_by, atoms_to_add), reverse=reverse)]
    a = molecule.GetAtom(oechem.OEHasMapIdx(sorted_atoms[0][0]))
    if 'ringsystem' in a.GetData():
        ring_idx = a.GetData('ringsystem')
        atoms.update(ring_systems[ring_idx][0])
        bonds.update(ring_systems[ring_idx][1])
    if 'fgroup' in a.GetData():
        fgroup = a.GetData('fgroup')
        atoms.update(functional_groups[fgroup][0])
        bonds.update(functional_groups[fgroup][1])
    atoms.add(sorted_atoms[0][0])
    bonds.add(sorted_atoms[0][1])
    return atoms, bonds

def get_atoms_bonds_from_fragment(fragments):
# Get atoms and bonds in map indices
    atoms_bonds = {}
    for bond_tuple in fragments:

        atoms = set()
        bonds = set()
        fragment_mol = fragmenter.fragments[bond_tuple]
        for atom in fragment_mol.GetAtoms():
            m = atom.GetMapIdx()
            if m < 1:
                continue
            atoms.add(m)
        for bond in fragment_mol.GetBonds():
            m1 = bond.GetBgn().GetMapIdx()
            m2 = bond.GetEnd().GetMapIdx()
            if m1 < 1 or m2 < 1:
                continue
            bonds.add((m1, m2))
        atoms_bonds[bond_tuple] = (atoms, bonds)
    return atoms_bonds

def get_wbo_growth(fragmenter, fragment_mol, atoms, bonds,  target_bond, heuristic):
    wbos = []
    fragments = []
    fragment_mol = fragmenter.fragments[target_bond]
    oechem.OEAddExplicitHydrogens(fragment_mol)
    try:
        charged_fragment = chemi.get_charges(fragment_mol, strict_stereo=False, strict_types=False)
        a1 = charged_fragment.GetAtom(oechem.OEHasMapIdx(target_bond[0]))
        a2 = charged_fragment.GetAtom(oechem.OEHasMapIdx(target_bond[-1]))
        bond = charged_fragment.GetBond(a1, a2)
        wbos.append(bond.GetData('WibergBondOrder'))
        fragments.append(charged_fragment)
    except RuntimeError:
        print('Bad fragment. Continuing to adding the next bond')
        pass
    while fragment_mol.GetMaxAtomIdx() < fragmenter.molecule.GetMaxAtomIdx():
        
        atoms, bonds = add_one_bond(molecule=fragmenter.molecule, atoms=atoms,
                                    bonds=bonds, functional_groups=fragmenter.functional_groups, 
                                    ring_systems=fragmenter.ring_systems, 
                                   target_bond=target_bond, heuristic=heuristic)

        ab_set = fragmenter._to_atom_bond_set(atoms, bonds)
        fragment_mol = fragmenter.atom_bond_set_to_mol(ab_set, expand_stereoisomers=False)
        oechem.OEAddExplicitHydrogens(fragment_mol)
        try:
            charged_fragment = chemi.get_charges(fragment_mol, strict_stereo=False, strict_types=False)
            a1 = charged_fragment.GetAtom(oechem.OEHasMapIdx(target_bond[0]))
            a2 = charged_fragment.GetAtom(oechem.OEHasMapIdx(target_bond[-1]))
            bond = charged_fragment.GetBond(a1, a2)
            fragment_wbo = bond.GetData('WibergBondOrder')
            wbos.append(fragment_wbo)
            mol_copy = copy.deepcopy(charged_fragment)
            fragments.append(mol_copy)
        except RuntimeError:
            pass
    return wbos, fragments


# Load in all kinase inhibitors
with open('../../combinatorial_fragmentation/filter/filtered_kinase_inhibitors.json', 'r') as f:
    kinase_inhibitors = json.load(f)
key_list = list(kinase_inhibitors.keys())

for ki in key_list:
    print(ki)
    try:
        os.mkdir(ki)
    except FileExistsError:
        print('{} directory already exists. Files will be overwritten'.format(ki))
    mapped_smiles = kinase_inhibitors[ki]['canonical_isomeric_explicit_hydrogen_mapped_smiles']
    chemi.mol_to_image_atoms_label(mapped_smiles, map_idx=True, fname='{}/{}_mapped.png'.format(ki, ki))
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, mapped_smiles)
    fragmenter = fragment.WBOFragmenter(mol)
    fragmenter.fragment(threshold=1.0, strict_stereo=False, strict_types=False, keep_non_rotor_ring_substituents=False) # Try without keeping nonrotors on ring
    atoms_bonds = get_atoms_bonds_from_fragment(fragmenter.fragments)
    fragments_to_drive = {}
    with PdfPages('{}/{}_fragment_growth_2.pdf'.format(ki, ki)) as pdf:
        for target_bond in fragmenter.fragments:
            key = workflow_api.serialize_key(target_bond)
            fragments_to_drive[key] = {'path_length': [], 'wbo': []}
            atoms = copy.deepcopy(atoms_bonds[target_bond][0])
            bonds = copy.deepcopy(atoms_bonds[target_bond][1])
            pl_wbos, pl_fragments = get_wbo_growth(fragmenter, fragmenter.fragments[target_bond], atoms, bonds,
                                                   target_bond, 'path_length')
            for fr, wb in zip(pl_fragments, pl_wbos):
                fragments_to_drive[key]['path_length'].append((oechem.OEMolToSmiles(fr), wb))
            pl_sizes = [f.GetMaxAtomIdx() for f in pl_fragments]
            atoms = copy.deepcopy(atoms_bonds[target_bond][0])
            bonds = copy.deepcopy(atoms_bonds[target_bond][1])
            wb_wbos, wb_fragments = get_wbo_growth(fragmenter, fragmenter.fragments[target_bond], atoms, bonds,
                                                  target_bond, 'wbo')
            wb_sizes = [f.GetMaxAtomIdx() for f in wb_fragments]
            for fr, wb in zip(wb_fragments, wb_wbos):
                fragments_to_drive[key]['wbo'].append((oechem.OEMolToSmiles(fr), wb))
            # plot

            plt.figure()
            plt.plot(pl_sizes, pl_wbos, color='salmon')
            plt.plot(pl_sizes, pl_wbos, 'o', color='salmon', label='path length')

            plt.plot(wb_sizes, wb_wbos, color='mediumslateblue')
            plt.plot(wb_sizes, wb_wbos, 'o', color='mediumslateblue', label='wbo')
            plt.legend()
            plt.xlabel('heavy atoms')
            plt.ylabel('Wiberg Bond Order')
            plt.title(target_bond)
            pdf.savefig()
            frags = copy.deepcopy(pl_fragments)
            for f in frags:
                for b in f.GetBonds():
                    b.DeleteData('WibergBondOrder')
            chemi.to_pdf(molecules=frags,bond_map_idx=target_bond, rows=3, cols=2, bo=pl_wbos, color=oechem.OELightSalmon,
                        oname='{}/{}_pl_frags_{}_{}_2.pdf'.format(ki, ki, target_bond[0], target_bond[1]))

            frags = copy.deepcopy(wb_fragments)
            for f in frags:
               for b in f.GetBonds():
                   b.DeleteData('WibergBondOrder')
            chemi.to_pdf(molecules=frags,bond_map_idx=target_bond, rows=3, cols=2, bo=wb_wbos, color=oechem.OELightBlue,
                       oname='{}/{}_wb_frags_{}_{}_2.pdf'.format(ki, ki, target_bond[0], target_bond[1]))

    with open('{}/{}_frags_to_drive_2.json'.format(ki, ki), 'w') as f:
       json.dump(fragments_to_drive, f, indent=2, sort_keys=True)






