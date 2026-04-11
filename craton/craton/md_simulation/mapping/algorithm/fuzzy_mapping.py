import functools
import itertools
import time
from collections import defaultdict
from typing import Dict, FrozenSet, Iterable, List, Optional, Set, Tuple

# other basic tools
import networkx
import numpy as np

from ....utils.geometry import calc_stru_para

# tools from compuchem

_MAX_DISTANCE = 2.4  # angstrom
_MIN_SIZE = 2


class TimeoutError(AssertionError):
    def __init__(self, msg, *args):
        self.default_result = (), ()
        super(TimeoutError, self).__init__(msg, *args)


def get_matches(
    st0,
    st1,
    max_distance=_MAX_DISTANCE,
    include_h=True,
    timeout=None,
    allowed_atoms0=None,
    allowed_atoms1=None,
    known_mapping=None,
    allow_ring_breaking=True,
    rotatable_bonds0=None,
    rotatable_bonds1=None,
    score_count_only=False,
):
    """
    本函数是查找MCS并获取MCS的主函数，其中的查找算法调用Mapper类中的方法
    Get a mapping between the atoms in st0 and st1 based on the positions of
    the atoms.

    Chirality, stereochemistry, and bond connectivity are also considered.

    :param st0: Structure to map against `st1`
    :type st0: schrodinger.structure.Structure
    :param st1: Structure to map against `st0`
    :type st1: schrodinger.structure.Structure

    :param max_distance: The maximum apart that heavy atoms will be in the
            mapping. Hydrogens may be farther than this value.
    :type distance: float (angstroms)

    :param include_h: Should hydrogen atoms be included in the map? If so,
            hydrogens are only allowed to map to other hydrogens.
    :type include_h: bool

    :param timeout: The max duration for which to attempt to find a mapping. If
            0, there is no timeout. If a timeout is hit, a TimeoutError is raised
    :type timeout: float (seconds)

    :param allowed_atoms0: What atoms in `st0` are allowed to participate in
            the map?
    :type allowed_atoms0: set of int (atom indices)
    :param allowed_atoms1: What atoms in `st1` are allowed to participate in
            the map?
    :type allowed_atoms1: set of int (atom indices)

    :param known_mapping: Atoms that are known to match. Keys are atom indices
            in `st0`, values are atom indices in `st1`
    :type known_mapping: dict {int: int}

    :param allow_ring_breaking: Are ring size changes and ring opening allowed?
    :type allow_ring_breaking: bool

    :param score_count_only: Score atom mapping only by atom counts?
    :type  score_count_only: bool

    :return: A dict of matched atom indices from st0 and st1
    :rtype: dict[int, int]

    :raise TimeoutError
    """
    rings0 = structure_find_rings(st0)
    rings1 = structure_find_rings(st1)

    # User passed in a known map
    if known_mapping and len(set(known_mapping.values())) != len(known_mapping):
        msg = "Input atom maps must include each atom only once. " "Duplicate atoms found for st1: {}".format(
            sorted(known_mapping.values())
        )
        raise ValueError(msg)

    mapper = Mapper(
        st0,
        st1,
        rings0,
        rings1,
        allowed_atoms0,
        allowed_atoms1,
        max_distance,
        False,
        allow_ring_breaking,
        known_mapping,
        rotatable_bonds0,
        rotatable_bonds1,
        score_count_only,
    )

    def mappings(timeout_point):
        """
        Build the map.

        If not initial guess, form a guess using only the atoms in the "core"
        of the structure, where core means heavy atoms.
        """
        yield from mapper.find_maps_iter(known_mapping=known_mapping, timeout_point=timeout_point)
        mapper.update_candidates(rings0, rings1, best_map, include_h)
        if include_h:  # Skip ring blob atoms' cis/trans check for Hydrogen atoms
            mapper.stereo_check_atoms0 = mapper.planar_atoms0
            mapper.stereo_check_atoms1 = mapper.planar_atoms1
        if best_map:
            yield from mapper.growMatchFromCore(best_map, timeout_point)
        else:  # in case core mapping fails (maybe it's unnecessary)
            yield from mapper.find_maps_iter(timeout_point=timeout_point)
        yield mapper.reassignTerminalAtoms(best_map)

    best_map = known_mapping or {}
    best_map_score = mapper.score(best_map)

    # Main mapping loop:
    timeout_point = timeout and time.time() + timeout
    for mapping in mappings(timeout_point):
        score = mapper.score(mapping)
        if score > best_map_score:
            best_map = mapping
            best_map_score = score
    
    best_map = ignore_match_ring_atom(st0,st1,best_map,known_mapping)
    return best_map

def ignore_match_ring_atom(st0,st1,mapd,known_mapping):
    remove_atoms = []
    reversed_mapd = {vv:kk for kk,vv in mapd.items()}
    matcheds_a = set(mapd.keys())
    matcheds_b = set(mapd.values())
    ring_blocks_a = {an:ii for ii,vv in enumerate(st0.ring_blocks) for an in vv}
    ring_blocks_b = {an:ii for ii,vv in enumerate(st1.ring_blocks) for an in vv}
    for ma,mb in mapd.items():
        if len(st0.Atoms[ma].has_ring) > 0 or len(st1.Atoms[mb].has_ring) > 0:
            if ma not in remove_atoms:
                atoma = st0.Atoms[ma]
                atomb = st1.Atoms[mb]
                hra = len(atoma.has_ring)
                hrb = len(atomb.has_ring)
                #ring_atoms_a = set([an  for rname in atoma.has_ring for an in st0._rings[rname][:-1]])
                #ring_atoms_b = set([an  for rname in atomb.has_ring for an in st1._rings[rname][:-1]])
                ring_atoms_a = st0.ring_blocks[ring_blocks_a[ma]] if ma in ring_blocks_a else []
                ring_atoms_a = set(ring_atoms_a)
                ring_atoms_b = st1.ring_blocks[ring_blocks_b[mb]] if mb in ring_blocks_b else []
                ring_atoms_b = set(ring_atoms_b)

                diffs_a = ring_atoms_a.difference(matcheds_a)
                diffs_b = ring_atoms_b.difference(matcheds_b)
                
                if len(diffs_a) > 0 or len(diffs_b) > 0 :
                    if len(diffs_a) > 0:
                        connects_atoms = [[an,ca] for ca in ring_atoms_a for an in st0.Atoms[ca].connectivity if an not in ring_atoms_a 
                                          ]
                        length_atoms = [st0.find_side_componend(*rn) + [rn[0]] for rn in connects_atoms]
                        match_length = [len(set(ln).intersection(matcheds_a)) for ln in length_atoms]
                        longest_index = match_length.index(max(match_length))
                        remove_atoms.extend(list(ring_atoms_a))
                        for ii, ans in enumerate(length_atoms):
                            if ii != longest_index:
                                remove_atoms.extend(ans)
                    else: #len(diffs_a) > 0:
                        connects_atoms = [[an,ca] for ca in ring_atoms_b for an in st1.Atoms[ca].connectivity if an not in ring_atoms_b 
                                          ]
                        length_atoms = [st1.find_side_componend(*rn) + [rn[0]] for rn in connects_atoms]
                        match_length = [len(set(ln).intersection(matcheds_b)) for ln in length_atoms]
                        longest_index = match_length.index(max(match_length))
                        reversed_atoms_a = [reversed_mapd[ann] for ann in ring_atoms_b if ann in reversed_mapd]
                        
                        remove_atoms.extend(reversed_atoms_a)
                        
                        for ii, ans in enumerate(length_atoms):
                            if ii != longest_index:
                                for anns in ans:
                                    if anns in reversed_mapd:
                                        remove_atoms.append(reversed_mapd[anns])
                        
    new_map = {}
    for ma,mb in mapd.items():
        if ma not in remove_atoms:
            new_map[ma] = mb
    return new_map

def old_ignore_match_ring_atom(st0,st1,mapd,known_mapping):
    remove_atoms = []
    matcheds_a = set(mapd.keys())
    matcheds_b = set(mapd.values())
    for ma,mb in mapd.items():
        if len(st0.Atoms[ma].has_ring) > 0:
            if ma not in remove_atoms:
                atoma = st0.Atoms[ma]
                atomb = st1.Atoms[mb]
                hra = len(atoma.has_ring)
                hrb = len(atomb.has_ring)
                ring_atoms_a = set([an  for rname in atoma.has_ring for an in st0._rings[rname][:-1]])
                ring_atoms_b = set([an  for rname in atomb.has_ring for an in st1._rings[rname][:-1]])
                diffs_a = ring_atoms_a.difference(matcheds_a)
                diffs_b = ring_atoms_b.difference(matcheds_b)
                if len(diffs_a) > 0 or len(diffs_b) > 0 :
                    connects_atoms = [[an,ca] for ca in ring_atoms_a for an in st0.Atoms[ca].connectivity if an not in ring_atoms_a 
                                      ]
                    length_atoms = [st0.find_side_componend(*rn) + [rn[0]] for rn in connects_atoms]
                    match_length = [len(set(ln).intersection(matcheds_a)) for ln in length_atoms]
                    longest_index = match_length.index(max(match_length))
                    remove_atoms.extend(list(ring_atoms_a))
                    for ii, ans in enumerate(length_atoms):
                        if ii != longest_index:
                            remove_atoms.extend(ans)
    new_map = {}
    for ma,mb in mapd.items():
        if ma not in remove_atoms:
            new_map[ma] = mb
    return new_map
        
def structure_as_graph(st):
    """
    Get a networkx Graph representing the basic topology of `st`. This does
    not include coordinates, bond orders, or atomic numbers. Nodes are atom
    indices. Terminal atoms (atoms with only 1 bond) are not included.

    :type st: schrodinger.structure.Structure
    :return: Graph of the input structure, nodes match input atom indices.
    :rtype: networkx.Graph
    """
    graph = networkx.Graph()
    for atom in st.Atoms:
        if len(atom.connect) == 1:
            continue
        atom_index = atom.No
        for neighbor in atom.connect:
            if len(st.Atoms[neighbor].connect) != 1:
                graph.add_edge(atom_index, st.Atoms[neighbor].No)
    return graph


class Mapper:
    """
    修改进度：
    Find a mapping between the atoms of st0 and st1.

    Mapping is primarily based on position, but it also takes into account
    stereochemistry and bonding.

    See find_maps_iter() and growMatchFromCore()

    各类mapping相关函数

    """

    def __init__(
        self,
        st0,
        st1,
        rings0=None,
        rings1=None,
        allowed_atoms0=None,
        allowed_atoms1=None,
        max_distance=_MAX_DISTANCE,
        include_hydrogens=False,
        allow_ring_breaking=True,
        known_mapping=None,
        rotatable_bonds0=None,
        rotatable_bonds1=None,
        score_count_only=False,
    ):
        if rings0 is None:
            rings0 = structure_find_rings(st0)
        if rings1 is None:
            rings1 = structure_find_rings(st1)

        if rotatable_bonds0 is None:
            self.rotatable_bonds0 = get_rotatable_bonds(st0, rings0)
        if rotatable_bonds1 is None:
            self.rotatable_bonds1 = get_rotatable_bonds(st1, rings1)

        # Bonds of the second structure
        self.st1_bonds = defaultdict(set)

        self.st0 = st0
        self.st1 = st1
        atom_groups0 = _get_atom_groups(st0, allowed_atoms0)
        atom_groups1 = _get_atom_groups(st1, allowed_atoms1)

        self.h_atoms0 = atom_groups0[1]
        self.h_atoms1 = atom_groups1[1]
        self.core_atoms0 = atom_groups0[0]
        self.core_atoms1 = atom_groups1[0]

        self.atom_matches = _get_atom_matches(
            st0,
            st1,
            rings0,
            rings1,
            atom_groups0,
            atom_groups1,
            max_distance=max_distance,
            include_hydrogens=include_hydrogens,
            allow_ring_breaking=allow_ring_breaking,
            known_mapping=known_mapping,
        )

        self.get_bond_matches(rings0, rings1, include_hydrogens)

        # Which atoms need to be checked for cis/trans?
        self.planar_atoms0 = _find_planar_atoms(st0)
        self.planar_atoms1 = _find_planar_atoms(st1)

        self.heavy_atoms0 = atom_groups0[0] | atom_groups0[2]
        self.heavy_atoms1 = atom_groups1[0] | atom_groups1[2]

        self.ring_atoms0 = set(itertools.chain.from_iterable(rings0))
        self.ring_atoms1 = set(itertools.chain.from_iterable(rings1))

        self._st0_dihedrals = {}
        self._st1_dihedrals = {}
        self.allow_ring_breaking = allow_ring_breaking

        self.blobs0, self.rings0 = _get_linked_ring_blobs(rings0)
        self.blobs1, self.rings1 = _get_linked_ring_blobs(rings1)

        ring_blob_atoms0 = frozenset().union(*self.blobs0)
        ring_blob_atoms1 = frozenset().union(*self.blobs1)

        # check cis/trans for ring blob atoms as well
        self.stereo_check_atoms0 = self.planar_atoms0.union(ring_blob_atoms0)
        self.stereo_check_atoms1 = self.planar_atoms1.union(ring_blob_atoms1)

        self.simple_rings0 = [r for r in self.rings0 if not r & ring_blob_atoms0]
        self.simple_rings1 = [r for r in self.rings1 if not r & ring_blob_atoms1]

        self.simple_ring_atoms0 = frozenset().union(*self.simple_rings0)
        self.simple_ring_atoms1 = frozenset().union(*self.simple_rings1)

        self.topo0 = structure_as_graph(st0)
        self.topo1 = structure_as_graph(st1)

        # FIXME: These `all_rings` include ring blobs with genus 2 as well.
        #        Not sure if they are still needed.
        self.all_rings0 = _get_all_ring_sets(rings0)
        self.all_rings1 = _get_all_ring_sets(rings1)

        self.known_mapping = known_mapping or {}
        self.score_count_only = score_count_only

    def get_bond_matches(self, rings0, rings1, include_hydrogens):

        bond_matches_raw = _get_bond_matches(self.st0, self.st1, self.atom_matches, rings0, rings1)
        # All allowed bond correspondences
        self.bond_matches = defaultdict(dict)
        for bond0, raw_bonds1 in bond_matches_raw.items():
            bonds1 = defaultdict(list)
            for bond in raw_bonds1:
                bonds1[bond[1]].append((bond[0], bond[2]))
                self.st1_bonds[bond[1]].add(bond[2])
            self.bond_matches[bond0[0]][bond0[1]] = bonds1

    def find_maps_iter(self, minimum_size=_MIN_SIZE, known_mapping=None, timeout_point=None):

        """
        Find all mappings between the two structures.

        @yield: dict of atom index in st0 -> atom index in st1
        """

        known_mapping = known_mapping or {}
        queue = []
        # Define starting point for the mapping in a queue
        for atom, candidates in self.atom_matches.items():
            if atom in known_mapping:
                continue
            candidates = tuple(((atom, c) for c in candidates))
            queue.append((known_mapping.copy(), candidates))
        return self._largerMaps_iter(queue, minimum_size=minimum_size, timeout_point=timeout_point)

    def growMatchFromCore(self, mapping, timeout_point=None):
        """
        Expand a known map, yielding each map that is larger.

        The atoms in the known map must be
        connected.

        @yield: dict of atom index in st0 -> atom index in st1
        """
        mapped = set(mapping.values())

        frontier = []
        for atom, st1_atom in mapping.items():
            for neighbor, allowed_bonds in self.bond_matches[atom].items():
                if neighbor in mapping:
                    continue
                candidates = [(s, neighbor, c) for (s, c) in allowed_bonds[st1_atom] if c not in mapped]
                frontier.extend(candidates)

        frontier.sort(reverse=True)
        frontier = tuple(((a, c) for (s, a, c) in frontier))
        queue = [(mapping.copy(), frontier)]
        return self._largerMaps_iter(queue, timeout_point=timeout_point)

    def _largerMaps_iter(self, queue, seen=None, minimum_size=_MIN_SIZE, timeout_point=None):
        """
        Generate maps that are larger than the input maps in queue. Does
        a breadth-first search to find the best mappings of atoms in
        self.st1 to the atoms in self.st2.

        The input queue is a list, where each element is:
        (mapping, frontier)

        The "frontier" is a list of:
        [st0_atom, st1_candidate_match]
        sorted from best match to worst match. The frontier can have multiple
        options for each atom, and those will all be attempted.

        """
        if seen is None:
            seen = set()

        def already_tried(mapping, frontier):
            """Don't retry something that we've already tried."""
            to_add = (frozenset(mapping.items()), frontier)
            if to_add not in seen:
                seen.add(to_add)
                return False
            return True

        while queue:
            mapping, frontier = queue.pop(0)
            back_mapping = {v: k for k, v in mapping.items()}
            while frontier:
                if timeout_point and time.time() > timeout_point:
                    raise TimeoutError("_largerMaps_iter")
                if already_tried(mapping, frontier):
                    mapping = tuple()
                    break

                next_frontier = []
                backtracking_was_added = False

                for i, (atom, candidate) in enumerate(frontier):
                    if atom in mapping or candidate in back_mapping:
                        continue
                    # if not self.bonding_ok(atom, candidate, mapping, back_mapping):
                    #     continue
                    if not self.stereo_ok(atom, candidate, mapping, back_mapping):
                        continue
                    if not backtracking_was_added:
                        backtracking_was_added = True
                        untried_frontier = frontier[i + 1 :] + ((atom, candidate),)
                        queue.append((mapping.copy(), untried_frontier))

                    mapping[atom] = candidate
                    back_mapping[candidate] = atom

                    for neighbor, allowed_bonds in self.bond_matches[atom].items():
                        if neighbor in mapping:
                            continue
                        for score, neighbor_candidate in allowed_bonds[candidate]:
                            if neighbor_candidate not in back_mapping:
                                next_frontier.append((score, neighbor, neighbor_candidate))

                next_frontier.sort(reverse=True)
                next_frontier = tuple(
                    ((a, c) for (s, a, c) in next_frontier if a not in mapping and c not in back_mapping)
                )
                frontier = next_frontier
            if len(mapping) > minimum_size:
                self.trimTerminalMismatches(mapping)
                if len(mapping) > minimum_size:
                    yield mapping

    def bothNeighbors(self, atom, candidate, mapping, back_mapping):
        """
        Return the union of two types of neighbors for `atom`
        1. Mapped neighbors in st0
        2. The corresponding atoms in st0 of mapped neighbors of `candidate` in st1
        :type         atom: int
        :type    candidate: int
        :type      mapping: dict of st0 atom indicies to st1 atom indices
        :type back_mapping: dict of st1 atom indicies to st0 atom indices

        :rtype: set of int
        """
        st0_nbs = {n for n in self.bond_matches[atom] if n in mapping}
        st1_nbs = (n for n in self.st1_bonds[candidate] if n in back_mapping)
        st0_nbs.update(map(back_mapping.get, st1_nbs))
        return st0_nbs

    def mappedNeighbors(self, atom, mapping):
        """The neighbors of an atom in st0"""
        return [n for n in self.bond_matches[atom] if n in mapping]

    def st1MappedNeighbors(self, atom, mapping):
        """The neighbors of an atom in st1"""
        return [n for n in self.st1_bonds[atom] if n in mapping]

    def stereo_ok(self, atom, candidate, current_map, back_map):
        """
        Would adding the mapping atom->candidate violate the stereochemistry
        of any neighbors of atom or candidate?

        :param        atom: Atom index of an atom in st0
        :type         atom: int
        :param   candidate: Atom index of an atom in st1 that may match `atom`
        :type    candidate: int
        :param current_map: Mapping of atom indices in st0 to st1
        :type  current_map: dict
        :param    back_map: Mapping of atom indices in st1 to st0
        :type     back_map: dict

        :return: Is the stereochemistry/chirality OK?
        :rtype: bool
        """
        both_neighbors = self.bothNeighbors(atom, candidate, current_map, back_map)
        for neighbor in both_neighbors:
            st1_neighbor = current_map[neighbor]
            next_neighbors = [n for n in self.bothNeighbors(neighbor, st1_neighbor, current_map, back_map) if n != atom]

            # Check the adjacent bond for cis/trans matching
            planar_neighbor = neighbor in self.stereo_check_atoms0
            planar_st1_neighbor = st1_neighbor in self.stereo_check_atoms1

            for next_neighbor in next_neighbors:
                st1_next_neighbor = current_map[next_neighbor]

                # Filter out next_neighbors that aren't equivalently bonded
                # to the candidate atom
                if st1_next_neighbor not in self.st1_bonds[st1_neighbor]:
                    continue

                # Only check for cis/trans across bonds where both atoms
                # are planar. This allows conversion from sp2 to sp3, while
                # preserving the orientation during that mutation.
                cis_trans_check = (planar_neighbor and next_neighbor in self.stereo_check_atoms0) or (
                    planar_st1_neighbor and st1_next_neighbor in self.stereo_check_atoms1
                )
                bond0 = tuple(sorted((neighbor, next_neighbor)))
                bond1 = tuple(sorted((st1_neighbor, st1_next_neighbor)))
                rot_bond_check = (
                    bond0 in self.rotatable_bonds0
                    and bond1 in self.rotatable_bonds1
                    and atom in self.core_atoms0
                    and candidate in self.core_atoms1
                )

                if not (cis_trans_check or rot_bond_check):
                    continue

                for a4 in self.bond_matches[next_neighbor]:
                    if a4 == neighbor:
                        continue
                    try:
                        st1_a4 = current_map[a4]
                    except KeyError:
                        continue
                    # Either both cis or both trans
                    st0_dih = self._dihedral(self.st0, atom, neighbor, next_neighbor, a4)
                    st1_dih = self._dihedral(self.st1, candidate, st1_neighbor, st1_next_neighbor, st1_a4)
                    abs_st0_dih = abs(st0_dih)
                    abs_st1_dih = abs(st1_dih)
                    if cis_trans_check:
                        if abs_st0_dih < 60 and abs_st1_dih > 120:
                            return False
                        if abs_st0_dih > 120 and abs_st1_dih < 60:
                            return False
                    if rot_bond_check:
                        mismatch = 260 > abs(st0_dih - st1_dih) > 100
                        if mismatch:
                            return False

            # Check the improper dihedrals for atoms with more than 2
            # connections.

            if (
                neighbor not in self.planar_atoms0
                and st1_neighbor not in self.planar_atoms1
                and len(next_neighbors) > 1
            ):

                next_neighbors = [neighbor] + next_neighbors
                for st0_atoms in itertools.combinations(next_neighbors, 3):
                    st1_atoms = list(map(current_map.get, st0_atoms))
                    st0_dih = self._dihedral(self.st0, atom, *st0_atoms)
                    st1_dih = self._dihedral(self.st1, candidate, *st1_atoms)
                    if st0_dih > 0 and st1_dih < 0 or st0_dih < 0 and st1_dih > 0:
                        return False
        return True

    def reassignTerminalAtoms(self, mapping, minimum_size=_MIN_SIZE, timeout_point=None):
        """
        Strips off the terminal atoms and adds them back.

        This ensures that all terminal atoms (atoms with one bond) are assigned
        to the best possible match. This is required because the mapping
        doesn't actually access every possible permutation.
        """
        atoms_with_terminal_neighbors = defaultdict(list)
        for atom_index in mapping:
            atom = self.st0.Atoms[atom_index]
            if len(atom.connect) == 1:
                parent = atom.connect[0]
                atoms_with_terminal_neighbors[parent].append(atom_index)
            elif atom_index not in atoms_with_terminal_neighbors:
                atoms_with_terminal_neighbors[atom_index] = []

        best_mapping = mapping
        best_score = self.score(mapping)
        seen = set()
        for parent, terminal_atoms in atoms_with_terminal_neighbors.items():
            # Remove these terminal atoms from the map
            reassigned = best_mapping.copy()
            for atom in terminal_atoms:
                del reassigned[atom]

            # Check all atoms attached to parent to see if they can be added.
            # _largerMaps_iter() will take care of checking the different
            # combinations.
            parent_match = mapping[parent]
            frontier = []
            for neighbor, allowed_bonds in self.bond_matches[parent].items():
                if neighbor in reassigned or len(self.st0.Atoms[neighbor].connect) != 1:
                    continue
                for score, neighbor_candidate in allowed_bonds[parent_match]:
                    frontier.append((score, neighbor, neighbor_candidate))
            frontier.sort(reverse=True)
            frontier = tuple(((a, c) for (s, a, c) in frontier))
            queue = [(reassigned, frontier)]

            for new_mapping in self._largerMaps_iter(queue, seen, minimum_size, timeout_point):
                new_score = self.score(new_mapping)
                if new_score > best_score:
                    best_score = new_score
                    best_mapping = new_mapping
        return best_mapping

    def _dihedral(self, st, a1, a2, a3, a4):
        """
        Get a dihedral from cache or structure. About a 20percent speedup.
        """

        if st == self.st0:
            cache = self._st0_dihedrals
        else:
            cache = self._st1_dihedrals

        try:
            return cache[(a1, a2, a3, a4)]
        except KeyError:
            value = calc_stru_para([st.Atoms[a1].coor, st.Atoms[a2].coor, st.Atoms[a3].coor, st.Atoms[a4].coor])
            cache[(a1, a2, a3, a4)] = value
            return value

    def _shouldRemoveTerminalAtom(self, atom_index0, mapping, topo0, topo1):

        """
        Should `atom_index0` be removed from `mapping`?

        Given atom_index0 which is terminal on structure 0, check its mapped
        atom is terminal on structure 1, and then remove the terminal atom from
        the map if one of the following is true:
        * The hybridization doesn't match
        * The element doesn't match
        * Either of the terminal atoms is a simple-ring atom (ring opening)

        """
        if atom_index0 in self.known_mapping:  # do not trim user input
            return False

        atom_index1 = mapping[atom_index0]
        # Only remove if the atom is terminal on both sides.
        neighbors1 = [n for n in self.st1MappedNeighbors(atom_index1, mapping.values()) if n in self.heavy_atoms1]
        if len(neighbors1) > 1:
            return False

        atom0 = self.st0.Atoms[atom_index0]
        atom1 = self.st1.Atoms[atom_index1]
        # element mismatch
        if atom0.atom_number != atom1.atom_number:
            return True
        # hybridization mismatch
        if len(atom0.connect) != len(atom1.connect):
            return True
        # partial terminal ring
        if _in_terminal_ring(
            atom_index0, self.simple_ring_atoms0, self.simple_rings0, topo0, mapping
        ) or _in_terminal_ring(atom_index1, self.simple_ring_atoms1, self.simple_rings1, topo1, mapping.values()):
            return True
        return False

    def removeTerminalAtoms(self, mapping, topo0, topo1):
        """
        See trimming rules in _shouldRemoveTerminalAtoms()
        """
        removed = set()
        outside0 = [n for n, d in topo0.degree() if d < 2]

        for atom_index0 in outside0:
            try:
                neighbors = list(topo0[atom_index0])
            except KeyError:
                # Already removed, or doesn't have any neighbors
                continue
            while len(neighbors) < 2:
                if not self._shouldRemoveTerminalAtom(atom_index0, mapping, topo0, topo1):
                    break
                removed.add(atom_index0)
                topo1.remove_node(mapping.pop(atom_index0))
                topo0.remove_node(atom_index0)
                if neighbors:
                    atom_index0 = neighbors[0]
                    neighbors = list(topo0[atom_index0])
                else:
                    break
        return removed

    def removeTerminalRings(self, mapping, ring_mapping, topo0, topo1):
        """
        Remove terminal ring if it's NOT mapped to an equal-sized ring.
        """
        to_remove = set()
        back_mapping = {v: k for k, v in mapping.items()}

        t_rings0 = _get_terminal_rings(topo0, self.rings0)
        t_rings1 = _get_terminal_rings(topo1, self.rings1)
        to_remove.update(itertools.chain.from_iterable(t_rings0 - set(ring_mapping)))
        to_remove.update(map(back_mapping.get, itertools.chain.from_iterable(t_rings1 - set(ring_mapping.values()))))
        for atom_index0 in to_remove:
            if atom_index0 not in self.known_mapping:
                topo0.remove_node(atom_index0)
                topo1.remove_node(mapping.pop(atom_index0))
        return to_remove

    def removeTerminalBlobs(self, mapping, ring_mapping, topo0, topo1):
        """
        remove if a linked ring system (blob) if both conditions are true
        1. it has only one external edge, i.e., terminal
        2. it has no substituent ring mapped to an equal-sized ring,
           i.e., terminal hopping
        """
        to_remove = set()
        back_mapping = {v: k for k, v in mapping.items()}
        atoms0, atoms1 = set(mapping), set(mapping.values())

        # candidate blobs: partially mapped and terminal
        c_blobs0 = {mapped for mapped in (b & atoms0 for b in self.blobs0) if mapped and _is_terminal(topo0, mapped)}
        c_blobs1 = {mapped for mapped in (b & atoms1 for b in self.blobs1) if mapped and _is_terminal(topo1, mapped)}
        to_remove.update(itertools.chain.from_iterable(b for b in c_blobs0 if all(not r <= b for r in ring_mapping)))
        to_remove.update(
            map(
                back_mapping.get,
                itertools.chain.from_iterable(b for b in c_blobs1 if all(not r <= b for r in ring_mapping.values())),
            )
        )
        for atom_index0 in to_remove:
            if atom_index0 not in self.known_mapping:
                topo1.remove_node(mapping.pop(atom_index0))
                topo0.remove_node(atom_index0)
        return to_remove

    def trimTerminalMismatches(self, mapping):
        """
        还未debug
        Remove mismatched terminal atoms, terminal rings, and terminal blobs.
        """
        # rings in the ring_mapping do not get trimmed
        ring_mapping = _get_ring_mapping(mapping, self.rings0, self.rings1)
        prev_count, removed = -1, set()  # 仅仅是设定一个终止条件

        topo0 = _subgraph_nx1(self.topo0, mapping)
        topo1 = _subgraph_nx1(self.topo1, mapping.values())

        while prev_count != len(removed):
            prev_count = len(removed)
            removed.update(self.removeTerminalAtoms(mapping, topo0, topo1))
            removed.update(self.removeTerminalRings(mapping, ring_mapping, topo0, topo1))
            removed.update(self.removeTerminalBlobs(mapping, ring_mapping, topo0, topo1))

        # Remove terminal atoms (H, F, Cl, etc) attached to removed atoms
        for atom in removed:
            for neighbor in self.bond_matches[atom]:
                if neighbor in self.topo0.nodes():
                    continue
                mapping.pop(neighbor, None)

    def update_candidates(
        self, rings0: Tuple[Tuple[int]], rings1: Tuple[Tuple[int]], best_map: Dict[int, int], include_h: bool
    ):
        """
        Update the atom candidates and bond candidates. Specifically
        * eliminate core (non-terminal) atom candidates using `best_map`
        * add ring atom candidates that are outside the distance cutoff
        * optionally add the H atom candidates
        * update bond matches
        """

        # use best_map to filter the atom_matches
        for index in self.atom_matches:
            if index in best_map:
                self.atom_matches[index] = [best_map[index]]

        # add out-of-range ring atom candidates
        ring_matches = _match_rings(self.rings0, self.rings1, best_map)
        for difference, r0, r1 in ring_matches:
            if difference:  # skip any mismatched ring pairs
                continue
            unmapped0 = r0 - best_map.keys()
            unmapped1 = r1 - set(best_map.values())
            for a0, a1 in itertools.product(unmapped0, unmapped1):
                if (
                    np.linalg.norm(
                        np.array(self.st0.Atoms[a0].coor)
                        - np.array(self.st1.Atoms[a1].coor)
                        # numpy.array(atom_xyz(self.st0.Atoms[a0])) - numpy.array(atom_xyz(self.st0.Atoms[a1]))
                    )
                    >= _MAX_DISTANCE
                ):
                    self.atom_matches[a0].append(a1)
        if include_h:  # add H atom candidates
            unmapped_h0 = self.h_atoms0 - set(best_map)
            unmapped_h1 = self.h_atoms1 - set(best_map.values())
            _update_hydrogen_matches(self.st0, self.st1, unmapped_h0, unmapped_h1, self.atom_matches)

        # TODO: this can be optimized by avoiding existing bond matches
        self.get_bond_matches(rings0, rings1, include_hydrogens=True)

        # update known_mapping to avoid unnecessary terminal trimming
        self.known_mapping = best_map

    def score(self, mapping):
        """
        Score a mapping based on RMSD and the number of atoms in the map.
        """

        if not mapping:
            return 0
        ALL_ATOM = 1
        HEAVY_ATOM = 1
        CHEMISTRY = 0.4
        ALL_RMSD = 0.2
        HEAVY_RMSD = 0.8
        RING_COUNT = 5
        RING_MATCH_NOT_RING = -1.5
        PARTIAL_RING_ATOMS = -0.1
        EXOCYCLIC_RING_ATOMS = -0.5
        RMSD_DEWEIGHT = 2.5

        atoms0, atoms1 = list(zip(*mapping.items()))

        heavy_atoms0 = [a for a in atoms0 if a in self.heavy_atoms0]
        heavy_atoms1 = [a for a in atoms1 if a in self.heavy_atoms1]

        # Atom count score
        # mapping原子数量占二者总原子数量的比例，用以表征两者相似程度
        all_atom_count_score = 2 * len(mapping) / (len(self.st0.Atoms) + len(self.st1.Atoms))

        # RMSD score
        all_rmsd = calc_in_place_rmsd(self.st0, list(atoms0), self.st1, list(atoms1))
        all_rmsd_score = RMSD_DEWEIGHT / (RMSD_DEWEIGHT + all_rmsd)

        if self.score_count_only:  # 预设的初始值为False，一般不用，但这里不知道什么时候会用到
            # use rmsd with small weight as tie breaker
            return all_atom_count_score + 0.01 * all_rmsd_score

        if self.heavy_atoms0 or self.heavy_atoms1:  # 用来比较0中heavy_atom是否比1中大？
            heavy_atom_count_score = 2 * len(heavy_atoms0) / float(len(self.heavy_atoms0) + len(self.heavy_atoms1))
        else:
            heavy_atom_count_score = 0

        if heavy_atoms0:  # 重原子的rmsd
            heavy_rmsd = calc_in_place_rmsd(self.st0, list(heavy_atoms0), self.st1, list(heavy_atoms1))
            heavy_rmsd_score = RMSD_DEWEIGHT / (RMSD_DEWEIGHT + heavy_rmsd)
        else:
            heavy_rmsd_score = 0

        # Chemistry score
        same_identity = 0
        for a0, a1 in mapping.items():
            a0 = self.st0.Atoms[a0]
            a1 = self.st1.Atoms[a1]
            if a0.atom_number == a1.atom_number:  # 检查原子序数是否相同，即元素种类
                same_identity += 1
            if len(a0.connect) == len(a1.connect):  # 该原子所连的键的数量是否相同
                same_identity += 1
        chemistry_score = same_identity / float(len(atoms0) + len(atoms1))  # 这块定义了化学相似性的打分用以判定

        # Ring similarity score
        # * If a ring is only partially mapped, we penalize using RMSD,(环被部分map下，不建议使用rmsd，为什么？)
        #   because this is often a ring flip (or at least it's worst when it
        #   is a ring flip)（环反转？）
        # * If an atom outside of a ring is inserted into a ring, the RMSD isn't
        #   what matters - it's basically just a matter of whether there is an
        #   outside atom inserted. (this is usually quite bad, and indicates
        #   a pretty obvious problem in the mapping).第二个原则也没太懂
        if self.ring_atoms0 or self.ring_atoms1:
            mapped1 = set(mapping.values())

            partial, exocyclic, matched_ring_count = _find_partial_ring_matches(
                mapping, self.all_rings0, self.all_rings1, self.ring_atoms0, self.ring_atoms1
            )
            ring_atom_count = float(len(self.ring_atoms0) + len(self.ring_atoms1))
            # More mapped rings is better!
            ring_count = matched_ring_count / ring_atom_count  # 对于多环分子，match上的环越多，match结果越好
            # In general, prefer mappings where rings atoms match to other ring atoms
            mismatched_ring_atoms = set([_f for _f in map(mapping.get, self.ring_atoms0) if _f is not None]) ^ (
                self.ring_atoms1 & mapped1
            )

            ring_match_not_ring = len(mismatched_ring_atoms) / ring_atom_count
            partial_ring_atoms = 0
            exocyclic_ring_atoms = 0
            if self.allow_ring_breaking:
                if partial:
                    partial = tuple(partial)
                    partial_ring_atoms = calc_in_place_rmsd(
                        self.st0, list(partial), self.st1, list(map(mapping.get, partial))
                    )
                if exocyclic:
                    exocyclic = tuple(exocyclic)
                    exocyclic_ring_atoms = calc_in_place_rmsd(
                        self.st0, list(exocyclic), self.st1, list(map(mapping.get, exocyclic))
                    ) * len(exocyclic)
            elif partial or exocyclic:
                partial_ring_atoms = 10000
                exocyclic_ring_atoms = 10000
        else:
            ring_count = 0
            ring_match_not_ring = 0
            partial_ring_atoms = 0
            exocyclic_ring_atoms = 0

        # RMSD is essentially a tie-breaker. We care much more about the size
        # of the map than about the atom positions. We care essentially not at
        # all about the position of hydrogens.
        return (
            all_atom_count_score * ALL_ATOM
            + heavy_atom_count_score * HEAVY_ATOM
            + chemistry_score * CHEMISTRY
            + ring_count * RING_COUNT
            + partial_ring_atoms * PARTIAL_RING_ATOMS
            + ring_match_not_ring * RING_MATCH_NOT_RING
            + exocyclic_ring_atoms * EXOCYCLIC_RING_ATOMS
            + all_rmsd_score * ALL_RMSD
            + heavy_rmsd_score * HEAVY_RMSD
        )

    # ==========================================================================
    def bonding_ok(self, atom, candidate, current_map, back_mapping):
        """
        If ring breaking is not allowed, check if both atoms have exactly the
        same neighbors.

        :param         atom: Atom index of an atom in st0
        :type          atom: int
        :param    candidate: Atom index of an atom in st1 that may match `atom`.
        :type     candidate: int
        :param  current_map: Mapping of atom indices in st0 to st1
        :type   current_map: dict
        :param back_mapping: Mapping of atom indices in st1 to st0
        :type  back_mapping: dict

        :return: Is the bonding/connectivity OK?
        :rtype: bool
        """
        if not self.allow_ring_breaking:
            existing_neighbors = self.mappedNeighbors(atom, current_map)
            candidate_neighbors = set(self.st1_bonds[candidate]) & set(back_mapping)
            atom_neighbors_mapped = set(map(current_map.get, existing_neighbors))
            return candidate_neighbors == atom_neighbors_mapped
        return True


def _get_atom_groups(
    ct, allowed_atoms: Optional[Set[int]] = None  # structure对象
) -> Tuple[Set[int], Set[int], Set[int], Set[int]]:
    """
    Return the core (not terminal), hydrogen, heavy terminal, and triple-bond
    atoms in the structure.

    对一个分子进行原子种类分组，返回四类原子
    传入：一个分子等参数
    返回：该分子中的四类原子
    """
    if allowed_atoms is None:
        allowed_atoms = set(range(0, len(ct.Atoms)))
    h_atoms = {a for a in allowed_atoms if ct.Atoms[a].atom_number == 1}
    core_atoms = allowed_atoms - h_atoms
    heavy_terminal_atoms = {a for a in core_atoms if len(ct.Atoms[a].connect) == 1}
    core_atoms -= heavy_terminal_atoms
    triple_bond_atoms = set()
    for atom in ct.Atoms:
        if "3" in atom.bond_type:
            triple_bond_atoms.add(atom.No)

    return core_atoms, h_atoms, heavy_terminal_atoms, triple_bond_atoms


def _get_atom_matches(
    ct0,
    ct1,
    rings0,
    rings1,
    atom_groups0,
    atom_groups1,
    max_distance=_MAX_DISTANCE,
    include_hydrogens=False,
    allow_ring_breaking=True,
    known_mapping=None,
):
    """
    * An atom match is allowed if either
        * it is in known_mapping
        * it is within max_distance angstroms.
    * If hydrogens are included in the match, hydrogens can match any other
      hydrogen bound to a heavy atom with the same element and number of bonds.

    :return: A dictionary of the allowed mappings for each atom
    :rtype: dict[int] -> list[int]
    """
    # use known_mapping directly as matches
    known_mapping = known_mapping or {}
    atom_matches = defaultdict(list, ((k, [v]) for k, v in known_mapping.items()))

    core_atoms0, h_atoms0, heavy_terminal_atoms0, triple_bond_atoms0 = atom_groups0
    core_atoms1, h_atoms1, heavy_terminal_atoms1, triple_bond_atoms1 = atom_groups1
    known0 = set(known_mapping)
    known1 = set(known_mapping.values())
    unmapped_core0 = core_atoms0 - known0
    unmapped_core1 = core_atoms1 - known1

    ring_atoms0 = set(itertools.chain.from_iterable(rings0))
    ring_atoms1 = set(itertools.chain.from_iterable(rings1))

    # Atom matches: spatially close and triple-bondness
    for atom_index in unmapped_core0:
        atom = ct0.Atoms[atom_index]

        matched_atoms = get_distance_ranking(atom, ct1, max_distance)
        matched_atoms = [m[0] for m in matched_atoms if m[0] in unmapped_core1]

        if triple_bond_atoms0 or triple_bond_atoms1:
            is_triple0 = atom_index in triple_bond_atoms0
            matched_atoms = [i for i in matched_atoms if (i in triple_bond_atoms1) == is_triple0]

        # Sort so that atoms that are in rings prefer to match to other
        # ring atoms, and non-ring atoms prefer to match non-ring atoms.
        # Ring/non-ring matches are allowed, they're just lower priority.
        #
        # Sort is stable, so this retains distance ordering within the
        # ring/non-ring groups.

        if atom_index in ring_atoms0:

            def key(a):
                return a in ring_atoms1

        else:

            def key(a):
                return a not in ring_atoms1

        if allow_ring_breaking:
            matched_atoms.sort(key=key)
        else:
            matched_atoms = list(filter(key, matched_atoms))

        if matched_atoms:
            atom_matches[atom.No] = matched_atoms

    # extend ring atom matches for small rings
    rings0 = [r for r in rings0 if len(r) < 7]
    rings1 = [set(r) for r in rings1 if len(r) < 7]
    if rings1:
        for r in rings0:
            unmatched0 = [i for i in r if i in core_atoms0 and not atom_matches[i]]
            candidates = set().union(*[set(atom_matches[i]) for i in r])
            overlaps = [candidates & r for r in rings1]
            max_overlap = max(len(o) for o in overlaps)
            if max_overlap == 0:
                break
            # The rings in st1 that overlap with r in st0 maximally are likely
            # the good ring candidates, thus extend the atom matches.
            for i, o in enumerate(overlaps):
                if len(o) != max_overlap:
                    continue
                unmatched1 = rings1[i] & core_atoms1 - o
                if not unmatched1 and not unmatched0:
                    continue
                for a0 in unmatched0:
                    for a1 in rings1[i]:
                        atom_matches[a0].append(a1)
                for a1 in unmatched1:
                    for a0 in r:
                        atom_matches[a0].append(a1)
    if include_hydrogens:
        _update_hydrogen_matches(ct0, ct1, h_atoms0 - known0, h_atoms1 - known1, atom_matches)
    atom_matches.update(_get_terminal_matches(ct0, ct1, heavy_terminal_atoms0 - known0, heavy_terminal_atoms1 - known1))
    return atom_matches


def _get_terminal_matches(st0, st1, terminal_atoms0, terminal_atoms1):
    """
    还没debug
    Find the atoms from `terminal_atoms1` in `st1` that could serve as
    matches for each atom in `terminal_atoms0` from `st0`.
    Here terminal atoms do not include hydrogens.

    :type st0: structure.Structure
    :type st1: structure.Structure
    :type terminal_atoms0: set of int
    :type terminal_atoms1: set of int
    :rtype: dict, int > [int]
    """
    # terminal_atoms0 = {1,2,3}
    # terminal_atoms1 = {4,5,6}
    atoms_by_number = defaultdict(list)
    for index in terminal_atoms1:
        atom = st1.Atoms[index]
        atoms_by_number[atom.atom_number].append(index)
    atom_matches = {}
    for index in terminal_atoms0:
        atom = st0.Atoms[index]
        candidates = atoms_by_number.get(
            atom.atom_number
        )  # terminal atom里面只有这一个准则，找原子序数相同的，也就是同种原子才可以进行terminal mapping

        if candidates:
            atom_matches[index] = candidates
    return atom_matches


def _get_bond_matches(ct0, ct1, atom_matches, rings0, rings1):
    """
    Get the allowed bond matches. Also deletes atoms from atom_matches if
    there is no bond match.

    returns a dict of:
    (st0_atom0, st0_atom1) -> list of (score, st1_atom0, st1_atom1)

    """
    allowed_atoms1 = set(itertools.chain.from_iterable(atom_matches.values()))

    # FIXME: Maybe bonds0 doesn't need both (i, j) and (j, i)

    bonds0 = _get_bond_vectors(ct0, atom_matches)
    bonds1 = _get_bond_vectors(ct1, allowed_atoms1)

    atoms0_with_bond_matches = set()

    # Bond matches: same orientation, starting from an allowed atom
    # sorted by best alignment
    bond_matches = {}
    for bond0 in bonds0:
        matches = []
        bond0_to_terminal_atom = len(ct0.Atoms[bond0[0][0]].connect) == 1 or len(ct0.Atoms[bond0[0][1]].connect) == 1
        bond0_in_ring = False if bond0_to_terminal_atom else _is_ring_bond(bond0[0], rings0)
        for bond1 in bonds1:
            if bond1[0][0] not in atom_matches[bond0[0][0]] or bond1[0][1] not in atom_matches[bond0[0][1]]:
                continue
            dotprod = np.dot(bond0[1], bond1[1])
            # Bonds are pointed in exact opposite directions
            if dotprod < -0.9 and not bond0_to_terminal_atom:
                continue
            bond1_in_ring = _is_ring_bond(bond1[0], rings1)

            score = dotprod
            if bond0_in_ring and bond1_in_ring:
                # Scoring only matters if there is a timeout or if we bail
                # after finding a certain number of bad matches.
                score += 0.5
            elif bond0_in_ring != bond1_in_ring:
                score -= 0.5
            matches.append((score, bond1[0][0], bond1[0][1]))
        if not matches:
            continue

        matches.sort(reverse=True)
        bond_matches[bond0[0]] = matches
        atoms0_with_bond_matches.update(bond0[0])

    extra = set(atom_matches) - atoms0_with_bond_matches
    if len(extra) > 1:  # to deal with edge cases like methane
        for atom in extra:
            del atom_matches[atom]
    return bond_matches


def _get_bond_vectors(ct, allowed_atoms: Iterable[int]) -> List[Tuple[Tuple[int, int], np.ndarray]]:
    """
    Find numpy array vectors along each bond.

    :rtype: list(tuple)
    :return: list of ((atom1, atom2), vector between them)
    """
    bonds = []
    for index in allowed_atoms:
        atom = ct.Atoms[index]
        for idx2 in atom.connect:
            atom2 = ct.Atoms[idx2]
            if atom2.No < atom.No:
                continue
            if atom2.No not in allowed_atoms:
                continue
            vector = _direction(atom, atom2)
            bonds.append(((atom.No, atom2.No), vector))
            bonds.append(((atom2.No, atom.No), -vector))
    return bonds


def _direction(atom0, atom1):
    vector = np.array(atom1.coor) - np.array(atom0.coor)
    vector /= np.linalg.norm(vector)
    return vector


def _is_ring_bond(bond, rings: Tuple[Tuple[int]]) -> bool:
    """
    还未debug
    :param bond:
    :param rings:
    :return:
    """
    # This ring bond detection depends on the fact that the ring atoms indices
    # in `rings` are ordered along bonds.
    a0, a1 = bond
    for ring in rings:
        if a0 in ring and a1 in ring:
            a0_i = ring.index(a0)
            a1_i = ring.index(a1)
            if abs(a0_i - a1_i) in (1, len(ring) - 1):
                return True
    return False


def _find_planar_atoms(st):
    """
    Find the planar atoms in `st`.
    planar atoms need to be checked for cis/trans stereochemistry.
    """
    planar_atoms = set()
    for atom in st.Atoms:
        bond_total = len(atom.connect)
        if bond_total == 2:
            for bond_type in atom.bond_type:
                if bond_type == "2":
                    planar_atoms.add(atom.No)
        elif bond_total == 3:
            n1, n2, n3 = atom.connect
            dih = calc_stru_para([st.Atoms[n1].coor, atom.coor, st.Atoms[n2].coor, st.Atoms[n3].coor])
            if abs(dih) > 150:
                planar_atoms.add(atom.No)
    return planar_atoms


def _get_linked_ring_blobs(rings: Tuple[Tuple[int]]) -> Tuple[Set[FrozenSet[int]], Set[FrozenSet[int]]]:
    """
    Merged rings with common nodes into blobs. Blobs include fused rings and
    spiro rings.
    """

    def get_edges(ring):
        ring_edges = list(zip(ring, ring[1:]))
        ring_edges.append((ring[-1], ring[0]))
        return ring_edges

    g = networkx.Graph()
    g.add_edges_from(itertools.chain.from_iterable(map(get_edges, rings)))
    ring_sets = {frozenset(r) for r in rings}
    blobs = {x for x in map(frozenset, networkx.connected_components(g)) if x not in ring_sets}
    return blobs, ring_sets


def _get_all_ring_sets(rings: Tuple[Tuple[int]]) -> Set[FrozenSet[int]]:
    """
    Given iterables of atom indices representing rings, find larger rings
    consisting of bridged ring systems.
    """
    all_rings = set()
    for r in rings:
        all_rings.add(frozenset(r))
    for r0, r1 in itertools.combinations(all_rings, 2):
        if len(r0 & r1) > 1:
            all_rings.add(r0 | r1)
    return all_rings


def _find_partial_ring_matches(
    mapping: Dict[int, int],
    rings0: Set[FrozenSet[int]],
    rings1: Set[FrozenSet[int]],
    ring_atoms0: Set[int],
    ring_atoms1: Set[int],
) -> Tuple[Set[int], Set[int], int]:
    """
    还未debug
    Find all atoms in mapping that are either:
    * in rings that are only partially mapped
    * in a ring in one structure, but not in a ring in the other structure.
    """
    back_map = {v: k for (k, v) in mapping.items()}
    mapped0 = set(mapping)
    mapped1 = set(mapping.values())

    ring_match_list = _match_rings(rings0, rings1, mapping)
    matched_ring_count = 0
    paired_ring_atoms = set()
    partial_ring_atoms = set()
    differing_atoms = set()
    for difference, r0, r1 in ring_match_list:
        if difference:
            # One or more atoms don't match the atoms in the other ring
            difference = [_f for _f in map(back_map.get, difference) if _f]
            differing_atoms.update(difference)
        if not r0.issubset(mapped0) and not r1.issubset(mapped1):
            partial_ring_atoms.update(r0 & mapped0)
        elif not difference:
            matched_ring_count += 1
            paired_ring_atoms.update(r0 & mapped0)

    # Do not consider ring atoms as exocyclic atoms
    to_exclude = {a for a in differing_atoms if a in ring_atoms0 and mapping[a] in ring_atoms1}

    exocyclic = differing_atoms - to_exclude
    return (partial_ring_atoms - paired_ring_atoms, exocyclic, matched_ring_count)


def _match_rings(rings0, rings1, mapping):
    """
    Figure out which rings in st1 correspond to which rings in st0. Returns
    full rings from each plus the atoms that don't match.

    :param  rings0: Sets of atom indices that represent rings
    :type   rings0: iterable of sets
    :param  rings1: Sets of atom indices that represent rings
    :type   rings1: iterable of sets
    :param mapping: Correspondence of atom indices in the first structure to
                    atom indices in the second structure.
    :type  mapping: dict
    :param  mapped: Mapped atoms in st1. Equivalent to set(mapping.values())
    :type   mapped: set

    :return: List of paired rings and the atoms that don't match. Format is:
           (mismatched atoms in rings1 indices, full ring from rings0, full ring
            from rings1) Indices in the difference correspond to st1.
    :rtype: List of tuples
    """
    mapped = frozenset(mapping.values())
    all_ring_matches = []
    for r0 in rings0:
        mapped_r0 = frozenset([_f for _f in map(mapping.get, r0) if _f])
        for r1 in rings1:
            if mapped_r0 & r1:
                difference = (mapped_r0 ^ r1) & mapped
                all_ring_matches.append((difference, r0, r1))

    # Sort the match list from most similar to least:
    all_ring_matches.sort(
        key=lambda d_r0_r1: (len(d_r0_r1[0]), len(d_r0_r1[0]) / float(len(d_r0_r1[1]) + len(d_r0_r1[2])))
    )

    # At this point all_ring_matches may contain multiple potential matches
    # for each ring. We want each ring to be uniquely assigned.
    matched0 = set()
    matched1 = set()
    # This will contain one and only one match for each ring
    best_ring_matches = []
    for d, r0, r1 in all_ring_matches:
        if not (r0 in matched0 or r1 in matched1):
            best_ring_matches.append((d, r0, r1))
            matched0.add(r0)
            matched1.add(r1)
    return best_ring_matches


def _update_hydrogen_matches(st0, st1, h_atoms0, h_atoms1, atom_matches):
    """
    还未debug
    For each hydrogen in st0, find the potential matching atoms in st1.

    Hydrogen atoms are allowed to match any hydrogen atom whose heavy atom is
    * allowed to match
    * of the same element and number of bonds

    :type atom_matches: `dict` of `int` to `list` of `int`
    """

    def get_attached_characteristics(st, h_atoms, allowed_atoms):
        attached_characteristics = {}
        host_to_H_map = defaultdict(list)
        for index in h_atoms:
            atom = st.Atoms[index]
            heavy = st.Atoms[atom.connect[0]]
            if heavy.No not in allowed_atoms:
                continue
            attached_characteristics[index] = ((len(heavy.connect), heavy.atom_number), heavy.No)
            host_to_H_map[heavy.No].append(index)
        return attached_characteristics, host_to_H_map

    allowed_atoms1 = set(itertools.chain.from_iterable(atom_matches.values()))
    st0_h, _ = get_attached_characteristics(st0, h_atoms0, atom_matches)
    st1_h, host_to_h1 = get_attached_characteristics(st1, h_atoms1, allowed_atoms1)
    for atom, characteristics in st0_h.items():
        host0 = characteristics[1]
        if host0 in atom_matches:
            hosts1 = atom_matches[host0]
            candidates = [
                h_index for host1 in hosts1 for h_index in host_to_h1[host1] if st1_h[h_index][0] == characteristics[0]
            ]
        if candidates:
            atom_matches[atom] = candidates


def _get_ring_mapping(mapping, rings0, rings1):
    """
    Get ring-to-ring match. Ring size change is NOT allowed.

    TODO: According to Yuqing Deng, we should also make sure the aromaticity
    match for sub-rings of linked rings.

    :param mapping: atom mapping
    :type  mapping: dict{int:int}
    :param  rings0: rings in molecule 0
    :type   rings0: set{frozenset{int}}
    :param  rings1: rings in molecule 1
    :type   rings1: set{frozenset{int}}

    :return: ring mapping
    :rtype : dict{frozenset{int}:frozenset{int}}
    """
    ring_mapping = {}
    for r0 in rings0:
        r1 = frozenset(map(mapping.get, r0))
        if r1 in rings1:
            ring_mapping[r0] = r1
    return ring_mapping


def _subgraph_nx1(g: networkx.Graph, nodes):
    """
    Extract subgraph.

    This is the equivalent implement in networkx 1. See DESMOND-9592
    """
    sg = g.__class__()
    sg.add_nodes_from(g.nbunch_iter(nodes))
    sg_adj = sg._adj
    g_adj = g._adj
    for n in sg.nodes:
        nbs = sg.adjlist_inner_dict_factory()
        sg_adj[n] = nbs
        for nb in g_adj[n]:
            if nb in sg_adj:
                sg_adj[nb][n] = None
                nbs[nb] = None
    return sg


def _in_terminal_ring(atom, ring_atoms, rings, topo, mapped):
    """
    还未debug
    :param       atom: atom index
    :param ring_atoms: simple ring atoms
    :type  ring_atoms: frozenset of int
    :param      rings: simple rings, i.e., rings not inside ring blobs
    :type       rings: list of frozenset of int
    :type        topo: networkx graph object
    :type      mapped: iterable
    """
    if atom not in ring_atoms:
        return False
    host_ring = [r for r in rings if atom in r][0]
    partial_ring = {a for a in host_ring if a in mapped}
    if _is_terminal(topo, partial_ring):
        return True


def _is_terminal(topo, mapped_blob):
    """
    还未debug
    A blob is terminal if its atoms have exactly 1 external edge.

    :param        topo: topology of the mapping
    :type         topo: networkx graph object
    :param mapped_blob: a set of mapped atoms
    :type  mapped_blob: frozenset{int}

    :rtype: bool
    """
    external_edge = 0
    for atom in mapped_blob:
        for neighbor in topo[atom]:
            if neighbor not in mapped_blob:
                external_edge += 1
        if external_edge > 1:
            return False
    else:
        return bool(external_edge)


def _get_terminal_rings(topo, rings):
    """
    还未debug
    Find terminal rings in a graph topology.

    :param  topo: topology of the mapping
    :type   topo: networkx graph object
    :param rings: rings in the molecule
    :type  rings: set{frozenset{int}}

    :return: terminal rings
    :rtype : set{frozenset{int}}
    """
    mapped_rings = [r for r in rings if r.issubset(topo)]
    t_filter = functools.partial(_is_terminal, topo)
    return set(filter(t_filter, mapped_rings))


def structure_find_rings(structure):
    return [ring[:-1] for ring in structure.ring_dict.values()]


def get_distance_ranking(atom0, st1, max_distance):
    atom_list = []
    for atom in st1.Atoms:
        distance = np.sum(np.square(np.array(atom0.coor) - np.array(atom.coor))) ** 0.5
        if distance > max_distance:
            continue
        atom_list.append((atom.No, distance))
    ranking_list = sorted(atom_list, key=lambda x: x[1])
    return ranking_list


def get_rotatable_bonds(molecule, rings, cutoff=6):
    from itertools import chain

    ring_atoms = set(chain.from_iterable(rings))
    total_heavy = sum(1 for a in molecule.Atoms if a.atom_number != 1)
    result = set()
    for bb in molecule.Bonds:
        if hasattr(bb, "flexible") and bb.flexible == "yes":
            if bb.a1 not in ring_atoms and bb.a2 not in ring_atoms:
                continue
            left_idx = molecule.find_side_componend(bb.a1, bb.a2)
            if len(left_idx) == len(molecule.Atoms) - 2:
                result.add((bb.a1, bb.a2))
            heavy_left = len([idx for idx in left_idx if molecule.Atoms[idx].atom_number != 1]) + 1
            heavy_right = total_heavy - heavy_left
            if (heavy_left > cutoff) and (heavy_right > cutoff):
                result.add((bb.a1, bb.a2))
    return result


def calc_in_place_rmsd(molecule0, atom_list0, molecule1, atom_list1):
    """
    calculate rmsd between two molecules
    :param molecule0: compuchem molecule object
    :param atom_list0: list of atom for molecule 0 that need to be involved in
    :param molecule1: compuchem molecule object
    :param atom_list1: list of atom for molecule 1 that need to be involved in
    :return: rmsd value for molecule 0 and 1
    """
    l = len(atom_list0)
    sum = 0
    for i in range(l):
        distance = 0
        for j in [0, 1, 2]:
            distance += (molecule0.Atoms[atom_list0[i]].coor[j] - molecule1.Atoms[atom_list1[i]].coor[j]) ** 2
        sum += distance

    rmsd = pow(sum / l, 0.5)
    return rmsd
