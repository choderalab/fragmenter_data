"""
Calculate improper angles of trivalent nitrogen if they are involved in a torsion scan.
"""

import json
import matplotlib.pyplot as plt
import numpy as np

import qcfractal.interface as ptl
from fragmenter import chemi
from openeye import oechem
import cmiles


def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)


def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


def find_improper_angles(mol):
    """
    Find the improper dihedral angles in some molecule. Currently supports
    those with a central trivalent nitrogen atom.
    Parameters
    ----------
    mol : OpenEye oemol
        oemol in which to look for improper angles
    Returns
    -------
    list
        Each element in the list is a 4-tuple of the coordinates for the
        atoms involved in the improper. The central atom is listed first
        in the tuple. Each member of the tuple is a numpy array.
    list
        List of strings for atoms in the improper, central atom is first.
    """
    mol_coords = mol.GetCoords()
    crdlist = []
    Idxlist = []
    for atom in mol.GetAtoms(oechem.OEIsInvertibleNitrogen()):
        # central atom
        aidx = atom.GetIdx()
        crd0 = np.asarray(mol_coords[aidx])
        # sort the neighbors
        nbors = sorted(list(atom.GetAtoms()))
        #check if there are 3 atoms connected to central atom in improper
        if len(nbors) != 3:
            return crdlist, namelist
        crd1 = np.asarray(mol_coords[nbors[0].GetIdx()])
        crd2 = np.asarray(mol_coords[nbors[1].GetIdx()])
        crd3 = np.asarray(mol_coords[nbors[2].GetIdx()])
        # store coordinates
        crdlist.append([crd0, crd1, crd2, crd3])
        Idxlist.append([atom.GetIdx(), nbors[0].GetIdx(), nbors[1].GetIdx(),nbors[2].GetIdx()])
    return crdlist, Idxlist


def calc_improper_angle(atom0, atom1, atom2, atom3, translate=True):
    """
    Calculate the improper dihedral angle of a set of given four atoms.
    Parameters
    ----------
    atom0 : numpy array
        CENTRAL atom coordinates
    atom1 : numpy array
        outer atom coordinates
    atom2 : numpy array
        outer atom coordinates
    atom3 : numpy array
        outer atom coordinates
    translate : bool
        True to translate central atom to origin, False to keep as is.
        This should not affect the results of the calculation.
    Returns
    -------
    float
        Angle in degrees.
    """
    if translate:
        atom1 = atom1 - atom0
        atom2 = atom2 - atom0
        atom3 = atom3 - atom0
        atom0 = atom0 - atom0 # central must be moved last
    # calculate vectors
    v0 = atom0-atom1
    v1 = atom2-atom1
    v2 = atom2-atom3
    w1 = np.cross(v0, v1)
    w2 = np.cross(v1, v2)
    angle = angle_between(w1,w2) # this angle should be in range [0,90]
    # compute distance from plane to central atom
    # eq 6 from http://mathworld.wolfram.com/Point-PlaneDistance.html
    # here I'm using atom1 for (x,y,z), but could also use atom2 or atom3
    numer = w2[0]*(atom0[0]-atom1[0]) + w2[1]*(atom0[1]-atom1[1]) + w2[2]*(atom0[2]-atom1[2])
    denom = np.sqrt(w2[0]**2 + w2[1]**2 + w2[2]**2)
    dist = numer/denom
    # set reference so that if central atom is above plane, angle -> [90,180]
    #print("this is the distance:" + str(dist))
    #if dist < 0:
    #    angle = angle*-1
    #angle = abs(angle)
    return angle

client = ptl.FractalClient()
ds = client.get_collection('TorsionDriveDataset', 'OpenFF Substituted Phenyl Set 1')

with open('data/qcarchive_torsiondrive_conformers.json', 'r') as f:
    conformers = json.load(f)
with open('data/qcarchive_torsiondrives.json', 'r') as f:
    energies = json.load(f)

for fgroup in conformers:
    print(fgroup)
    confs = conformers[fgroup]
    energy = energies[fgroup]
    all_angles = []
    for i, index in enumerate(energy['indices']):
        angles = []
        oemols = [cmiles.utils.mol_from_json(c) for c in confs[energy['indices'][i]]]
        fname = 'data/{}_conformers_{}.xyz'.format(fgroup, i)
        ofs = oechem.oemolostream(fname)
        ofs.SetFormat(oechem.OEFormat_XYZ)
        for mol in oemols:
            oechem.OEWriteMolecule(ofs, mol)

        entry = ds.get_entry(index)
        dih = entry.td_keywords.dihedrals[0]
        # First, figure out the index needed
        crds, indxs = find_improper_angles(oemols[0])
        if not crds:
            print('no trivalent nitrogen in {}'.format(fgroup))
            continue
        idx_to_use = None
        for j, idx in enumerate(indxs):
            n = 0
            for k in idx:
                if k in dih:
                    n += 1
            if n == 3:
                idx_to_use = j
                print(idx, dih)
        if idx_to_use == None:
            print('trivalent nitrogen is not in dihedral {}'.format(fgroup))
            continue
        for mol in oemols:
            crds, idx = find_improper_angles(mol)
            crds = crds[idx_to_use]
            angle = calc_improper_angle(crds[0], crds[1], crds[2], crds[3])
            angles.append(angle * 180 / np.pi)
        all_angles.append(angles)

    plt.figure()
    colors = chemi._KELLYS_COLORS
    angles = np.arange(-165, 195, 15)
    for i, angle in enumerate(all_angles):
        plt.plot(angles, angle, 'o', color=colors[i], label=i)
        plt.plot(angles, angle, color=colors[i], label=i, linewidth=0.8)
        plt.title('Improper angles, {}'.format(fgroup), fontsize=14)
        plt.xlabel('Torsion angle (degrees)', fontsize=14)
        plt.ylabel('Improper angles (degrees)', fontsize=14)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14);
        plt.savefig('figures/qcarchive_torsiondrives/{}_improper_angles.pdf'.format(fgroup))