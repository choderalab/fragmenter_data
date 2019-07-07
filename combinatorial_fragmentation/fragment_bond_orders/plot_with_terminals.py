import fragmenter
import json
import cmiles
import qcfractal.interface as portal
import oenotebook as oenb
import glob
from openeye import oechem, oequacpac
import seaborn as sbn
import joypy
import matplotlib.pyplot as plt
import os
from matplotlib.backends.backend_pdf import PdfPages


with open('../filter/filtered_kinase_inhibitors.json', 'r') as f:
    kinase_inhibitors = json.load(f)
# Remove Trametinib because it didn't complete yet
kinase_inhibitors.pop('Trametinib')
kinase_inhibitors.pop('Abemaciclib')

for name in kinase_inhibitors:
    print(name)
    with open('{}/{}_oe_wbo.json'.format(name, name), 'r') as f:
        oe_wbo = json.load(f)
    with open('{}/{}_parent_oe_wbo.json'.format(name, name), 'r') as f:
        oe_wbo_parent = json.load(f)
    #with open('{}/{}_oe_wbo_by_bond.json'.format(name, name), 'r') as f:
    #    frag_with_bond = json.load(f)
    with open('{}/{}_oe_wbo_by_bond_with_parent.json'.format(name, name), 'r') as f:
        frag_with_bond = json.load(f)

    # deserialize
    oe_wbo_deserialized = {}
    for frag in oe_wbo:
        oe_wbo_deserialized[frag] = {}
        for bond in oe_wbo[frag]:
            key = (bond.split('[')[-1].split(']')[0].split(','))
            key = (int(key[0]), int(key[-1]))
            oe_wbo_deserialized[frag][key] = oe_wbo[frag][bond]

    oe_wbo_parent_deserialized = {}
    for bond in oe_wbo_parent:
        key = key = (bond.split('[')[-1].split(']')[0].split(','))
        key = (int(key[0]), int(key[-1]))
        oe_wbo_parent_deserialized[key] = oe_wbo_parent[bond]

    # combine with rest of fragments
    oe_wbo_deserialized['parent'] = oe_wbo_parent_deserialized

    # load in fragments
    with open('../fragment/kinase_inhitibors/{}_fragments.json'.format(name), 'r') as f:
        fragments = json.load(f)

    frag_with_bond = {}
    for frag in oe_wbo_deserialized:
        if frag == 'parent':
            map_to_parent = kinase_inhibitors[name]['canonical_isomeric_explicit_hydrogen_mapped_smiles']
        else:
            map_to_parent = fragments[frag]['provenance']['routine']['enumerate_fragments']['map_to_parent']
        oemol = cmiles.utils.load_molecule(map_to_parent)
        r_bonds = []
        for bond in oemol.GetBonds():
            a1 = bond.GetBgn()
            a2 = bond.GetEnd()
            if a1.IsHydrogen() or a2.IsHydrogen():
                continue
            if a1.IsHalogen() or a2.IsHalogen():
                continue
            if bond.GetOrder() > 1:
                continue
            if bond.IsInRing():
                continue
            key = (a1.GetMapIdx(), a2.GetMapIdx())
            if key not in frag_with_bond and tuple(reversed(key)) not in frag_with_bond:
                frag_with_bond[key] = {}
            r_bonds.append(key)
        #r_bonds = [(b.GetBgn().GetMapIdx(), b.GetEnd().GetMapIdx()) for b in oemol.GetBonds() if b.IsRotor()]
        for b in r_bonds:
            terminal = False
            # check if it's terminal
            bond = oemol.GetBond(oemol.GetAtom(oechem.OEHasMapIdx(b[0])), oemol.GetAtom(oechem.OEHasMapIdx(b[1])))
            if not bond.IsRotor():
                # This is terminal
                terminal = True
            reverse = tuple(reversed(b))
            if b not in oe_wbo_deserialized[frag]:
                # Try reverse
                bo = oe_wbo_deserialized[frag][reverse]
            else:
                bo = oe_wbo_deserialized[frag][b]
            if terminal:
                frag_key = frag + '_term'
            else:
                frag_key = frag
            try:
                frag_with_bond[b][frag_key] = bo
            except KeyError:
                frag_with_bond[reverse][frag_key] = bo

    # sort fragments by OE WBO
    sorted_frags = {}
    for b in frag_with_bond:
        list_1 = []
        list_2 = []
        for frag in frag_with_bond[b]:
            list_1.append(frag)
            list_2.append(frag_with_bond[b][frag]['ensamble'])
            sorted_frags[b] =  [x for _,x in sorted(zip(list_2,list_1))]
    rot_bonds = list(frag_with_bond.keys())

    with PdfPages('{}/{}_fragment_bond_order_with_terminal.pdf'.format(name, name)) as pdf:
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
                elif frag.split('_')[-1] == 'term':
                    sbn.kdeplot(wbo, shade= True, color='green', alpha=0.8)
                else:
                    sbn.kdeplot(wbo, shade= True, alpha=0.8)
                sbn.distplot(wbo, hist=False, rug=True, kde=False, color='black')
                sbn.kdeplot(wbo, lw=1, color='black')
                plt.axvline(x=wbo_s, ymin=0, ymax=1, color='black', linewidth=0.5)
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

    # hihglight bond in fragments
    for b in sorted_frags:
        dir_name = '{}/bond_{}_{}'.format(name, b[0], b[1])
        try:
            os.mkdir(dir_name)
        except FileExistsError:
            pass
        for i, frag in enumerate(sorted_frags[b]):
            terminal = False
            if frag.split('_')[-1] == 'term':
                frag = frag.split('_')[0]
                terminal = True
            try:
                bo = oe_wbo_deserialized[frag][b]['ensamble']
            except KeyError:
                bo = oe_wbo_deserialized[frag][(b[1], b[0])]['ensamble']
            if frag == 'parent':
                map_from_parent = kinase_inhibitors[name]['canonical_isomeric_explicit_hydrogen_mapped_smiles']
                fname = '{}/parent_{}_bond_{}_{}.png'.format(dir_name, i, b[0], b[-1])
            else:
                map_from_parent = fragments[frag]['provenance']['routine']['enumerate_fragments']['map_to_parent']
            #b_1 = b.split('[')[-1].split(']')[0].split(',')
            #b_2 = (int(b_1[0]), int(b_1[1]))
            if terminal:
                fname = '{}/frag_{}_term_bond_{}_{}.png'.format(dir_name, i, b[0], b[1])
            elif not frag == 'parent':
                fname = '{}/frag_{}_bond_{}_{}.png'.format(dir_name, i, b[0], b[1])
            fragmenter.chemi.highlight_bond_by_map_idx(map_from_parent, [b], wbo=bo, fname=fname)



