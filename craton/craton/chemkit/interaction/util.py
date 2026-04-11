import itertools
from collections import namedtuple

import numpy as np
##from openbabel import pybel

from .config import *


def ring_is_planar(r_atoms):
    """Given a set of ring atoms, check if the ring is sufficiently planar
    to be considered aromatic"""
    normals = []
    n = len(r_atoms)
    for i in range(n):
        # Check for neighboring atoms in the ring
        left_atom = (i - 1) % n
        right_atom = (i + 1) % n
        coor_this = r_atoms[i].coor
        coor_left = r_atoms[left_atom].coor
        coor_right = r_atoms[right_atom].coor
        normals.append(np.cross((coor_left - coor_this), (coor_right - coor_this)))
    # Given all normals of ring atoms and their neighbors, the angle between any has to be 5.0 deg or less
    for n1, n2 in itertools.product(normals, repeat=2):
        arom_angle = vecangle(n1, n2)
        if all([arom_angle > AROMATIC_PLANARITY, arom_angle < 180.0 - AROMATIC_PLANARITY]):
            return False
    return True


def ions_not_in_ring(r_atoms):
    for atom in r_atoms:
        if atom.atom_name.upper() in METAL_IONS:
            return False
    return True


def vecangle(v1, v2, deg=True):
    if np.array_equal(v1, v2):
        return 0.0
    dm = np.dot(v1, v2)
    cm = np.linalg.norm(v1) * np.linalg.norm(v2)
    angle = np.arccos(dm / cm)  # Round here to prevent floating point errors
    return np.degrees([angle])[0] if deg else angle


def centroid(coo) -> np.ndarray:
    return np.array(list(map(np.mean, (([c[0] for c in coo]), ([c[1] for c in coo]), ([c[2] for c in coo])))))


def normalize_vector(v):
    norm = np.linalg.norm(v)
    return v / norm if not norm == 0 else v


def projection(pnormal1, ppoint, tpoint):
    """Calculates the centroid from a 3D point cloud and returns the coordinates
    :param pnormal1: normal of plane
    :param ppoint: coordinates of point in the plane
    :param tpoint: coordinates of point to be projected
    :returns : coordinates of point orthogonally projected on the plane
    """
    # Choose the plane normal pointing to the point to be projected
    pnormal2 = np.array([coo * (-1) for coo in pnormal1])
    d1 = euclidean3d(tpoint, pnormal1 + ppoint)
    d2 = euclidean3d(tpoint, pnormal2 + ppoint)
    if d1 is not None and d2 is not None:
        pnormal = pnormal1 if d1 < d2 else pnormal2
    else:
        return None
    # Calculate the projection of tpoint to the plane
    sn = -np.dot(pnormal, vector(ppoint, tpoint))
    sd = np.dot(pnormal, pnormal)
    sb = sn / sd
    return [c1 + c2 for c1, c2 in zip(tpoint, [sb * pn for pn in pnormal])]


def euclidean3d(v1, v2):
    """Faster implementation of euclidean distance for the 3D case."""
    return np.sqrt((v1[0] - v2[0]) ** 2 + (v1[1] - v2[1]) ** 2 + (v1[2] - v2[2]) ** 2)


def vector(p1, p2):
    """Vector from p1 to p2.
    :param p1: coordinates of point p1
    :param p2: coordinates of point p2
    :returns : numpy array with vector coordinates
    """
    return None if len(p1) != len(p2) else np.array([p2[i] - p1[i] for i in range(len(p1))])


def cluster_doubles(double_list):
    """Given a list of doubles, they are clustered if they share one element
    :param double_list: list of doubles
    :returns : list of clusters (tuples)
    """
    location = {}  # hashtable of which cluster each element is in
    clusters = []
    # Go through each double
    for t in double_list:
        a, b = t[0], t[1]
        # If they both are already in different clusters, merge the clusters
        if a in location and b in location:
            if location[a] != location[b]:
                if location[a] < location[b]:
                    clusters[location[a]] = clusters[location[a]].union(clusters[location[b]])  # Merge clusters
                    clusters = clusters[: location[b]] + clusters[location[b] + 1 :]
                else:
                    clusters[location[b]] = clusters[location[b]].union(clusters[location[a]])  # Merge clusters
                    clusters = clusters[: location[a]] + clusters[location[a] + 1 :]
                # Rebuild index of locations for each element as they have changed now
                location = {}
                for i, cluster in enumerate(clusters):
                    for c in cluster:
                        location[c] = i
        else:
            # If a is already in a cluster, add b to that cluster
            if a in location:
                clusters[location[a]].add(b)
                location[b] = location[a]
            # If b is already in a cluster, add a to that cluster
            if b in location:
                clusters[location[b]].add(a)
                location[a] = location[b]
            # If neither a nor b is in any cluster, create a new one with a and b
            if not (b in location and a in location):
                clusters.append(set(t))
                location[a] = len(clusters) - 1
                location[b] = len(clusters) - 1
    return map(tuple, clusters)


def is_lig(hetid):
    """Checks if a PDB compound can be excluded as a small molecule ligand"""
    h = hetid.upper()
    return not (h == "HOH" or h in UNSUPPORTED)


def nucleotide_linkage(residues):
    """Support for DNA/RNA ligands by finding missing covalent linkages to stitch DNA/RNA together."""

    nuc_covalent = []
    #######################################
    # Basic support for RNA/DNA as ligand #
    #######################################
    nucleotides = ["A", "C", "T", "G", "U", "DA", "DC", "DT", "DG", "DU"]
    dna_rna = {}  # Dictionary of DNA/RNA residues by chain
    covlinkage = namedtuple("covlinkage", "id1 chain1 pos1 conf1 id2 chain2 pos2 conf2")
    # Create missing covlinkage entries for DNA/RNA
    for ligand in residues:
        resname, chain, pos = ligand
        if resname in nucleotides:
            if chain not in dna_rna:
                dna_rna[chain] = [
                    (resname, pos),
                ]
            else:
                dna_rna[chain].append((resname, pos))
    for chain in dna_rna:
        nuc_list = dna_rna[chain]
        for i, nucleotide in enumerate(nuc_list):
            if not i == len(nuc_list) - 1:
                name, pos = nucleotide
                nextnucleotide = nuc_list[i + 1]
                nextname, nextpos = nextnucleotide
                newlink = covlinkage(
                    id1=name,
                    chain1=chain,
                    pos1=pos,
                    conf1="",
                    id2=nextname,
                    chain2=chain,
                    pos2=nextpos,
                    conf2="",
                )
                nuc_covalent.append(newlink)

    return nuc_covalent


def int32_to_negative(int32):
    """Checks if a suspicious number (e.g. ligand position) is in fact a negative number represented as a
    32 bit integer and returns the actual number.
    """
    dct = {}
    if int32 == 4294967295:  # Special case in some structures (note, this is just a workaround)
        return -1
    for i in range(-1000, -1):
        dct[np.uint32(i)] = i
    if int32 in dct:
        return dct[int32]
    else:
        return int32


def sort_members_by_importance(members):
    """Sort the members of a composite ligand according to two criteria:
    1. Split up in main and ion group. Ion groups are located behind the main group.
    2. Within each group, sort by chain and position."""
    main = [x for x in members if x[0] not in METAL_IONS]
    ion = [x for x in members if x[0] in METAL_IONS]
    sorted_main = sorted(main, key=lambda x: (x[1], x[2]))
    sorted_main = sorted(main, key=lambda x: (x[1], x[2]))
    sorted_ion = sorted(ion, key=lambda x: (x[1], x[2]))
    return sorted_main + sorted_ion


def classify_by_name(names):
    """Classify a (composite) ligand by the HETID(s)"""
    if len(names) > 3:  # Polymer
        if len(set(RNA).intersection(set(names))) != 0:
            ligtype = "RNA"
        elif len(set(DNA).intersection(set(names))) != 0:
            ligtype = "DNA"
        else:
            ligtype = "POLYMER"
    else:
        ligtype = "SMALLMOLECULE"

    for name in names:
        if name in METAL_IONS:
            if len(names) == 1:
                ligtype = "ION"
            else:
                if "ION" not in ligtype:
                    ligtype += "+ION"
    return ligtype


def canonicalize(lig, preserve_bond_order=False):
    """Get the canonical atom order for the ligand."""
    atomorder = None
    # Get canonical atom order

    lig = pybel.ob.OBMol(lig.OBMol)
    if not preserve_bond_order:
        for bond in pybel.ob.OBMolBondIter(lig):
            if bond.GetBondOrder() != 1:
                bond.SetBondOrder(1)
    lig.DeleteData(pybel.ob.StereoData)
    lig = pybel.Molecule(lig)
    testcan = lig.write(format="can")
    try:
        pybel.readstring("can", testcan)
        reference = pybel.readstring("can", testcan)
    except IOError:
        return atomorder
    if testcan != "":
        reference.removeh()
        isomorphs = get_isomorphisms(reference, lig)  # isomorphs now holds all isomorphisms within the molecule
        if not len(isomorphs) == 0:
            smi_dict = {}
            smi_to_can = isomorphs[0]
            for x in smi_to_can:
                smi_dict[int(x[1]) + 1] = int(x[0]) + 1
            atomorder = [smi_dict[x + 1] for x in range(len(lig.atoms))]
        else:
            atomorder = None
    return atomorder


def get_isomorphisms(reference, lig):
    """Get all isomorphisms of the ligand."""
    query = pybel.ob.CompileMoleculeQuery(reference.OBMol)
    mappr = pybel.ob.OBIsomorphismMapper.GetInstance(query)
    if all:
        isomorphs = pybel.ob.vvpairUIntUInt()
        mappr.MapAll(lig.OBMol, isomorphs)
    else:
        isomorphs = pybel.ob.vpairUIntUInt()
        mappr.MapFirst(lig.OBMol, isomorphs)
        isomorphs = [isomorphs]
    return isomorphs
