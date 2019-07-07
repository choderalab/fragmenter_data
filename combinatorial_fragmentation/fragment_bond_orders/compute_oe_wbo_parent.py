import fragmenter
import json
import cmiles
from openeye import oechem, oequacpac
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import seaborn as sbn


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="calculate OE WBO")
    parser.add_argument('-n', '--name', type=str, help='Molecule name')
    args = parser.parse_args()
    name = args.name

    with open('../filter/filtered_kinase_inhibitors.json', 'r') as f:
        kinase_inhibitors = json.load(f)
    kinase_inhibitors[name] = kinase_inhibitors[name]

    mapped_smiles = kinase_inhibitors[name]['canonical_isomeric_explicit_hydrogen_mapped_smiles']
    oemol = cmiles.utils.load_molecule(mapped_smiles, toolkit='openeye')
    charged = fragmenter.chemi.get_charges(oemol, keep_confs=-1)

    oe_wbo_full = {}
    for bond in charged.GetBonds():
        bond_key = (bond.GetBgn().GetMapIdx(), bond.GetEnd().GetMapIdx())
        oe_wbo_full[bond_key] = {'ensamble': bond.GetData('WibergBondOrder')}
        oe_wbo_full[bond_key]['individual_conf'] = []

    for i, conf in enumerate(charged.GetConfs()):
        mol_copy = oechem.OEMol(conf)
        # Get WBO
        if oequacpac.OEAssignPartialCharges(mol_copy, oequacpac.OECharges_AM1BCCSym):
            for bond in mol_copy.GetBonds():
                bond_key = (bond.GetBgn().GetMapIdx(), bond.GetEnd().GetMapIdx())
                try:
                    oe_wbo_full[bond_key]['individual_conf'].append(bond.GetData('WibergBondOrder'))
                except KeyError:
                    if 0 in bond_key:
                        oe_wbo_full[bond_key] = {'individual_conf': [bond.GetData('WibergBondOrder')]}
                    else:
                        reverse_key = tuple(reversed(bond_key))
                        oe_wbo_full[reverse_key]['individual_conf'].append(bond.GetData('WibergBondOrder'))
        else:
            print('AM1BCC charging failed for {}, {}'.format(str(i), i))

    # serialize and save
    serialized = {}
    for bond in oe_wbo_full:
        key = fragmenter.workflow_api.serialize_key(bond)
        serialized[key] = oe_wbo_full[bond]
    # save file
    with open('{}/{}_parent_oe_wbo.json'.format(name, name), 'w') as f:
        json.dump(serialized, f, indent=2, sort_keys=True)

    # replot to include a red distribution for parent
    # load others, deserialize and replot
    with open('{}/{}_oe_wbo_by_bond.json'.format(name, name), 'r') as f:
        by_bond = json.load(f)
    frag_with_bond = {}
    for bond in by_bond:
        key = bond.split('[')[-1].split(']')[0].split(',')
        key = (int(key[0]), int(key[1]))
        frag_with_bond[key] = by_bond[bond]
    # add parent
    for bond in frag_with_bond:
        try:
            full = oe_wbo_full[bond]
        except KeyError:
            key = (bond[-1], bond[0])
            full = oe_wbo_full[key]
        frag_with_bond[bond]['parent'] = full
    # serialize and save

    serialized_with_parent = {}
    for bond in frag_with_bond:
        key =fragmenter. workflow_api.serialize_key(bond)
        serialized_with_parent[key] = frag_with_bond[bond]
    with open('{}/{}_oe_wbo_by_bond_with_parent.json'.format(name, name), 'w') as f:
        json.dump(serialized_with_parent, f, indent=2, sort_keys=True)

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
    with PdfPages('{}/{}_fragment_bond_orders_with_parent.pdf'.format(name, name)) as pdf:
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
                if frag == 'parent':
                    sbn.kdeplot(wbo, shade= True, color='red', alpha=0.8)
                else:
                    sbn.kdeplot(wbo, shade= True,  alpha=0.8)
                sbn.distplot(wbo, hist=False, rug=True, kde=False, color='black')
                sbn.kdeplot(wbo, lw=1, color='black')
                plt.axvline(x=wbo_s, ymin=0, ymax=1, color='black', linewidth=0.5)

                plt.xlim(x_min-0.05, x_max+0.05)
                plt.yticks([])
                ax.yaxis.set_label_coords(-0.05, 0)
                plt.ylabel(i, rotation=0, size=8)
                if i != n-1:
                    plt.xticks([])
                else:
                    plt.xlabel('Bond order')
                if i == 0:
                    plt.title('bond {}'.format(b))
            overlap=0.5
            h_pad = 5 + (- 5*(1 + overlap))
            fig.tight_layout(h_pad=h_pad)
            pdf.savefig(bbox_inches='tight')
            plt.close()

    for b in frag_with_bond:
        try:
            wbo = oe_wbo_full[b]['ensamble']
        except KeyError:
            wbo = oe_wbo_full[(b[-1], b[0])]['ensamble']
        fragmenter.chemi.highlight_bond_by_map_idx(mapped_smiles, [b], wbo=wbo, fname='{}/parent_bond_{}_{}.png'.format(name, b[0], b[1]))




