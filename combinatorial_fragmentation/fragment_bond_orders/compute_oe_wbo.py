import fragmenter
import json
import cmiles
from openeye import oechem, oequacpac
import os
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import seaborn as sbn


def compute_oe_wbo(molecule_input):
    oe_wbo = {}
    for frag in molecule_input:
        print(frag)
        oe_wbo[frag] = {}
        # ensamble WBO
        map_to_parent = molecule_input[frag]['provenance']['routine']['enumerate_fragments']['map_to_parent']
        oemol = cmiles.utils.load_molecule(map_to_parent)
        charged = fragmenter.chemi.get_charges(oemol)
        for bond in charged.GetBonds():
            bond_key = (bond.GetBgn().GetMapIdx(), bond.GetEnd().GetMapIdx())
            oe_wbo[frag][bond_key] = {'ensamble': bond.GetData('WibergBondOrder')}
            oe_wbo[frag][bond_key]['individual_conf'] = []
        for i, json_mol in enumerate(molecule_input[frag]['initial_molecule']):
            oemol = cmiles.utils.mol_from_json(json_mol)
            # Add map to parent to keep track of bonds
            #map_to_parent = ruxolitinib_inputs[frag]['provenance']['routine']['enumerate_fragments']['map_to_parent']
            atom_map_to_parent = cmiles.utils.get_atom_map(oemol, map_to_parent, strict=False)
            invert_map = cmiles.utils.invert_atom_map(atom_map_to_parent)
            # Add parent map to frag
            for atom in oemol.GetAtoms():
                idx = atom.GetIdx()
                if idx in invert_map:
                    atom.SetMapIdx(invert_map[atom.GetIdx()])
            #AM1 calculations
            oemol_copy = oechem.OEMol(oemol)
            #am1charges = oequacpac.OEAM1Charges(optimize=False)
            if oequacpac.OEAssignPartialCharges(oemol_copy, oequacpac.OECharges_AM1BCCSym):
                # Save WBO
                for bond in oemol_copy.GetBonds():
                    bond_key = (bond.GetBgn().GetMapIdx(), bond.GetEnd().GetMapIdx())
                    try:
                        oe_wbo[frag][bond_key]['individual_conf'].append(bond.GetData('WibergBondOrder'))
                    except KeyError:
                        if 0 in bond_key:
                            oe_wbo[frag][bond_key] = {'individual_conf': [bond.GetData('WibergBondOrder')]}
                        else:
                            reverse_key = tuple(reversed(bond_key))
                            oe_wbo[frag][reverse_key]['individual_conf'].append(bond.GetData('WibergBondOrder'))
            else:
                print('AM1BCC charging failed for {}, {}'.format(str(i), frag))
    return oe_wbo


def serialize(oe_wbo):
    serialized_wbo = {}
    for frag in oe_wbo:
        serialized_wbo[frag] = {}
        for bond in oe_wbo[frag]:
            key = fragmenter.workflow_api.serialize_key(bond)
            serialized_wbo[frag][key] = oe_wbo[frag][bond]
    return serialized_wbo

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="calculate OE WBO")
    parser.add_argument('-i', '--infile', type=str, help='Input JSON conformers file')
    parser.add_argument('-n', '--name', type=str, help='Molecule name')
    args = parser.parse_args()
    infile = args.infile
    name = args.name

    with open('../filter/filtered_kinase_inhibitors.json', 'r') as f:
        kinase_inhibitors = json.load(f)
    kinase_inhibitors[name] = kinase_inhibitors[name]

    with open(infile, 'r') as f:
        fragments_inputs = json.load(f)
    oe_wbo = compute_oe_wbo(fragments_inputs)
    serialized = serialize(oe_wbo)
    try:
        os.mkdir(name)
    except FileExistsError:
        print('{} directory already exists. Files will be overwritten'.format(name))
    with open('{}/{}_oe_wbo.json'.format(name, name), 'w') as f:
        json.dump(serialized, f, indent=2, sort_keys=True)

    # organize for easier plotting
    # organize so it's easier to make plots
    mapped_parent_smiles = kinase_inhibitors[name]['canonical_isomeric_explicit_hydrogen_mapped_smiles']
    mapped_parent_mol = cmiles.utils.load_molecule(mapped_parent_smiles)
    fragmenter.chemi.mol_to_image_atoms_label(mapped_parent_smiles, '{}/mapped_{}.png'.format(name, name))
    rot_bonds = []
    for bond in mapped_parent_mol.GetBonds():
        if bond.IsRotor():
            key = (bond.GetBgn().GetMapIdx(), bond.GetEnd().GetMapIdx())
            rot_bonds.append(key)
    frag_with_bond = {b:{} for b in rot_bonds }
    for frag in oe_wbo:
        map_to_parent = fragments_inputs[frag]['provenance']['routine']['enumerate_fragments']['map_to_parent']
        oemol = cmiles.utils.load_molecule(map_to_parent)
        r_bonds = [(b.GetBgn().GetMapIdx(), b.GetEnd().GetMapIdx()) for b in oemol.GetBonds() if b.IsRotor()]
        for b in r_bonds:
            reverse = tuple(reversed(b))
            if b not in oe_wbo[frag]:
                # Try reverse
                bo = oe_wbo[frag][reverse]
            else:
                bo = oe_wbo[frag][b]
            try:
                frag_with_bond[b][frag] = bo
            except KeyError:
                frag_with_bond[reverse][frag] = bo
    serialized = {}
    for bond in frag_with_bond:
        key = fragmenter.workflow_api.serialize_key(bond)
        serialized[key] = frag_with_bond[bond]
    with open('{}/{}_oe_wbo_by_bond.json'.format(name, name), 'w') as f:
        json.dump(serialized, f, indent=2, sort_keys=True)

    # sort fragments by wbo
    sorted_frags = {}
    for b in frag_with_bond:
        list_1 = []
        list_2 = []
        for frag in frag_with_bond[b]:
            list_1.append(frag)
            list_2.append(frag_with_bond[b][frag]['ensamble'])
            sorted_frags[b] = [x for _,x in sorted(zip(list_2, list_1))]

    rot_bonds = list(frag_with_bond.keys())

    # plot results on one pdf page
    with PdfPages('{}/{}_fragment_bond_orders.pdf'.format(name, name)) as pdf:
        for b in rot_bonds:
            #b = rot_bonds[3]
            n = len(frag_with_bond[b])

            fig, axes = plt.subplots(n, 1)
            fig.dpi = 400
            x_min = 3
            x_max = 0
            for f in frag_with_bond[b]:
                wbo = frag_with_bond[b][f]['individual_conf']
                if min(wbo) < x_min:
                    x_min = min(wbo)
                if max(wbo) > x_max:
                    x_max = max(wbo)

            for i, frag in enumerate(sorted_frags[b]):
                wbo = frag_with_bond[b][frag]['individual_conf']

                wbo_s = frag_with_bond[b][frag]['ensamble']
                ax = plt.subplot(n, 1, i+1)
                ax.spines['left'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.spines['top'].set_visible(False)
                ax.spines['bottom'].set_visible(False)
                ax.patch.set_facecolor('none')
                sbn.kdeplot(wbo, shade= True, alpha=0.8)
                sbn.distplot(wbo, hist=False, rug=True, kde=False, color='black')
                sbn.kdeplot(wbo, lw=1, color='black')
                plt.axvline(x=wbo_s, ymin=0, ymax=1, color='black', linewidth=0.5)
                #if i == 0:
                    #sbn.kdeplot(wbo)
                    #sbn.kdeplot(bond_order_std_aromatic_unopt['Imatinib']['wiberg_bo'][bond], shade=True, label='non-optimized')
                    #try:
                    #    key = (bond[0]-1, bond[1]-1)
                    #    sbn.kdeplot(oe_wbo[key], shade=True, label='openeye')
                    #except KeyError:
                    #    key = (bond[-1]-1, bond[0]-1)
                    #    sbn.kdeplot(oe_wbo[key], shade=True, label='openeye')
                #else:
                #    sbn.kdeplot(bond_order_std_aromatic_opt['Imatinib']['wiberg_bo'][bond], shade=True)
                #    sbn.kdeplot(bond_order_std_aromatic_unopt['Imatinib']['wiberg_bo'][bond], shade=True)
                #    try:
                 #       key = (bond[0]-1, bond[1]-1)
                 #       sbn.kdeplot(oe_wbo[key], shade=True)
                  #  except KeyError:
                  #      key = (bond[-1]-1, bond[0]-1)
                  #      sbn.kdeplot(oe_wbo[key], shade=True)
                plt.xlim(x_min-0.05, x_max+0.05)
                plt.yticks([])
                #ax.set_yticklabels(bond_order_std_rot_bonds_opt['Imatinib']['bonds'], fontsize=6, rotation=0)
                #ax.yaxis.grid(False)
                ax.yaxis.set_label_coords(-0.05, 0)
                plt.ylabel(i, rotation=0, size=8)
                if i != n-1:
                    plt.xticks([])
                else:
                    plt.xlabel('Bond order')
                if i == 0:
                    #plt.legend(prop={'size': 10}, bbox_to_anchor=(1.35, 1))
                    plt.title('bond {}'.format(b))
            overlap=0.5
            h_pad = 5 + (- 5*(1 + overlap))
            fig.tight_layout(h_pad=h_pad)
            pdf.savefig(bbox_inches='tight')
            plt.close()

    for b in sorted_frags:
        for i, frag in enumerate(sorted_frags[b]):
            try:
                bo = oe_wbo[frag][b]['ensamble']
            except KeyError:
                bo = oe_wbo[frag][(b[1], b[0])]['ensamble']
            map_from_parent = fragments_inputs[frag]['provenance']['routine']['enumerate_fragments']['map_to_parent']
            #b_1 = b.split('[')[-1].split(']')[0].split(',')
            #b_2 = (int(b_1[0]), int(b_1[1]))
            fragmenter.chemi.highlight_bond_by_map_idx(map_from_parent, [b], wbo=bo,
                                                       fname='{}/frag_{}_bond_{}_{}.png'.format(name, i, b[0], b[1]))




