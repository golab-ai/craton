import networkx as nx
import numpy


class FragmentGroup:
    """
    Groups of mutating fragments, that attach to at least one common core atom
    """

    def __init__(self):
        # list of sets of atoms in reactant
        # the list index is the fragment id in this group
        self._source_atoms = []
        # dictionary of attachment points keyed by fragment index
        self._source_att_points = {}
        # list of sets of atoms in reactant
        # the list index is the fragment id in this group
        self._dest_atoms = []
        # dictionary of attachment points keyed by fragment index
        self._dest_att_points = {}

        # list of attachment point tuples ((sa, sb), (da, db))
        # where sa is the core atom in reactant;
        # where da is the core atom in product.
        # (sa, da) is a matching atom pair
        # sb is a dummy atom connected to sa in reactant
        # db is a dummy atom connected to da in product
        # both sb and db can be None
        self._attachment_points = []

    def addAttachmentPoint(self, attachment_point):
        """
        add one pair of attachment points ((sa, sb), (da, db))
        attachment_point: two-element tuple for both source and dest structures,
                          (core_atom_idx, bridgeatom_idx)
        """
        if attachment_point not in self._attachment_points:
            self._attachment_points.append(attachment_point)

    def addAttachmentPoints(self, attachment_points):
        """
        add a bunch of attachmentpoints
        attachment_points: list of attachment_point
        """
        for ap in attachment_points:
            self.addAttachmentPoint(ap)

    def addSourceAtoms(self, source_atoms):
        """
        add source fragment
        source_atoms: set of source atom indices
        """
        if not frozenset(source_atoms) in self._source_atoms:
            self._source_atoms.append(frozenset(source_atoms))

    def addDestAtoms(self, dest_atoms):
        """
        add dest fragment
        dest_atoms: set of dest atom indices
        """
        if not frozenset(dest_atoms) in self._dest_atoms:
            self._dest_atoms.append(frozenset(dest_atoms))

    def addAttpointToSourceFragment(self, s_frag_id, att_point):
        """
        populate the attachment point dictionary for source fragment,
                 s_frag_id
        s_frag_id: index of fragment in the source fragment list
        att_point: attachment point in the form of (sa, sb)
                   sa is the core atom where the fragment is connected;
                   sb is the dummy atom in the fragment that is connected
                      to sa
        """
        if s_frag_id in self._source_att_points:
            self._source_att_points[s_frag_id].add(att_point)
        else:
            self._source_att_points[s_frag_id] = set()
            self._source_att_points[s_frag_id].add(att_point)

    def addAttpointToDestFragment(self, d_frag_id, att_point):
        """
        populate the attachment point dictionary for dest fragment,
                 d_frag_id
        d_frag_id: index of fragment in the source fragment list
        att_point: attachment point in the form of (da, db)
                   da is the core atom where the fragment is connected;
                   db is the dummy atom in the fragment that is connected
                      to da
        """
        if d_frag_id in self._dest_att_points:
            self._dest_att_points[d_frag_id].add(att_point)
        else:
            self._dest_att_points[d_frag_id] = set()
            self._dest_att_points[d_frag_id].add(att_point)

    def getSourceFragmentAttpoints(self, s_frag_id):
        """
        get attachment points for source fragment, s_frag_id
        s_frag_id: index into source fragment list
        return: list of attachment points
        """
        return self._source_att_points.get(s_frag_id, [])

    def getDestFragmentAttpoints(self, d_frag_id):
        """
        get attachment points for dest fragment, d_frag_id
        d_frag_id: index into dest fragment list
        return: list of attachment points
        """
        return self._dest_att_points.get(d_frag_id, [])

    def getAttachmentPoints(self):
        """
        get all attachmentpoints as a list of tuple, ((sa, sb), (da, db))
        """
        return self._attachment_points

    def getSourceAtoms(self):
        """
        return list of source fragment atoms
        """
        return self._source_atoms

    def getNumSourceFrags(self):
        """
        return: number of fragments in reactant
        """
        return len(self._source_atoms)

    def getDestAtoms(self):
        """
        return list of dest fragment atoms
        """
        return self._dest_atoms

    def getNumDestFrags(self):
        """
        return: number of fragments in product
        """
        return len(self._dest_atoms)

    def getSourceAttachmentPoints(self):
        """
        get attachment points for all source fragments
        """
        ret_val = []
        for ((sa, sb), (da, db)) in self._attachment_points:
            if sb:
                ret_val.append((sa, sb))
        return ret_val

    def getDestAttachmentPoints(self):
        """
        get attachment points for all dest fragments
        """
        ret_val = []
        for ((sa, sb), (da, db)) in self._attachment_points:
            if db:
                ret_val.append((da, db))
        return ret_val


def find_connnected_atoms(atom_idx, fragment, ct):
    """
    find all atom a in fragment that is connected to atom_idx
    atom_idx: atom index
    fragment: list of atom indices
    ct: structure
    return list of atom indices  in fragment that is connected to atom_idx
    """

    ret_val = []
    atom = ct.Atoms[atom_idx]
    for idx in fragment:
        frag_atom = ct.Atoms[idx]
        if atom.No in frag_atom.connect:
            ret_val.append(idx)

    return ret_val


def convert_fragment_to_graph(bonds, atoms):
    """
    Convert a molecule to networkx graph, with graph properties
    of atoms that should not be included in soft bond
    :param bonds: all bonds in the molecule with bond order
    :type bonds: Dict[(Int, Int):Int]
    :param atoms: list of fragment atoms
    :type atoms: List[Int]
    :rtype: networkx.Graph
    """
    g = nx.Graph(dont_cut=set())
    g.add_nodes_from(atoms)

    for (a, b), v in bonds.items():
        if (a in atoms) and (b in atoms):
            g.add_edge(a, b, order=v)
            if v == 3:
                g.graph["dont_cut"].add(a)
                g.graph["dont_cut"].add(b)

    return g


def group_connected_atoms(atoms, ct):
    """
    group connected atom in to lists of lists
    atoms: atom indices
    ct: structure
    return: list of sets, each atom list is a connected fragment
    """

    # import networkx as nx

    # bonds = {(x.atom1.index, x.atom2.index): x.order for x in ct.bond}
    # g = convert_fragment_to_graph(bonds, atoms)
    g = nx.Graph()
    g.add_nodes_from(atoms)
    for atom in atoms:
        for neighbor in ct.Atoms[atom].connect:
            if neighbor in atoms:
                g.add_edge(atom, neighbor)
    return [frozenset(x) for x in nx.connected_components(g)]


def find_connected_fragments(atom_index, fragment_list, ct):
    """
    find all connected fragments to given atom_index
    """
    conn_fragments = []
    for (frag_idx, frag) in enumerate(fragment_list):
        conn_atoms = find_connnected_atoms(atom_index, frag, ct)
        if conn_atoms:
            conn_fragments.append((frag_idx, conn_atoms))
    return conn_fragments


def generate_connection_points(source_att, source_bridge_list, dest_att, dest_bridge_list):
    """
    source_att: source core atom index
    source_bridge_list: list of atoms in source fragment that is attached to source core atom
    dest_att: dest core atom index
    dest_bridge_list: list of atoms in dest fragment that is attached to source core atom
    return: list of attachment points ((sa, ab), (da, db) ...)
    """

    min_length = min(len(source_bridge_list), len(dest_bridge_list))

    att_points = []
    for i in range(min_length):
        att_points.append(((source_att, source_bridge_list[i]), (dest_att, dest_bridge_list[i])))

    for i in range(min_length, len(source_bridge_list)):
        att_points.append(((source_att, source_bridge_list[i]), (dest_att, None)))
    for i in range(min_length, len(dest_bridge_list)):
        att_points.append(((source_att, None), (dest_att, dest_bridge_list[i])))

    return att_points


def register_att_point_fragmentgroup(att_point_in_frag_group, att, conn_atoms, group_id, frag_id):
    """
    store the group_id and frag_id for attachment point
    att_point_in_frag_group: dictionary keyed by att_point atom pair, input and output
    att: core atom attachment point
    conn_atoms: atom is fragment bonded to att
    group_id: index to fragmentgroup
    group_id: index to fragment within a fragment group
    """

    for f_atom in conn_atoms:
        if (att, f_atom) not in att_point_in_frag_group:
            att_point_in_frag_group[(att, f_atom)] = (group_id, frag_id)


def update_visted_fragments(frag_visited, frag_id, other_frag_id, anchor_atom, connected_atoms):
    """
    update visted fragments, add other_id to list of other fragments, add attachment points
    frag_visited: input and output dictionary for visited fragments
    frag_id: id of fragment
    other_frag_id: id of fragment sharing anchor atom
    anchor_atom: anchor atom index
    connected_atoms: all atoms connected to the anchor atom
    return: Nothing
    """

    if str(frag_id) in frag_visited:
        if str(anchor_atom) in frag_visited[str(frag_id)][1]:
            frag_visited[str(frag_id)][1][str(anchor_atom)].extend(connected_atoms)
        else:
            frag_visited[str(frag_id)][1][str(anchor_atom)] = connected_atoms
    else:
        frag_visited[str(frag_id)] = [[], {str(anchor_atom): connected_atoms}]
    if other_frag_id != -1:
        frag_visited[str(frag_id)][0].append(other_frag_id)


def find_all_related_fragments(source_frag_id, source_frag_visited, dest_frag_visited):
    """
    find all fragments in both source and dest cts that share attachment point
    source_frag_id: source fragment input - hash index into source_frag_visited
    source_frag_visited: source fragment attachment information
    dest_frag_visited: dest fragment attachment information
    return: [sets of source fragment ids, sets of dest fragment ids]
    """

    ret_val = [set([str(source_frag_id)]), set()]

    for d_frag in source_frag_visited[source_frag_id][0]:
        ret_val[1].add(str(d_frag))
        for s_frag in dest_frag_visited[str(d_frag)][0]:
            ret_val[0].add(str(s_frag))

    return ret_val


def get_mutating_fragment_groups(source_core_atoms, dest_core_atoms, source_ct, dest_ct):
    """
    given matched source and dest core atoms, cts, find all fragment groups sharing
        attachment anchor atoms in the core.  Two dictionaries keyed by
        source and dest attachment points are also returned to locate attachmen point
        in which fragment of which fragment group
    source_core_atoms: lists of core atoms for initial molecule
    dest_core_atoms: lists of core atoms for initial molecule
    source_ct: structure of initial molecule
    dest_ct: structure of final molecule
    return: (lists FragmentGroup objects,
             dictionary of tuples (group number, source fragment number)
                        keyed by source attachment points,
             dictionary of tuples (group number, dest fragment number)
                        keyed by source attachment points,
    """

    source_non_core = set(range(len(source_ct.Atoms))) - set(source_core_atoms)
    dest_non_core = set(range(len(dest_ct.Atoms))) - set(dest_core_atoms)

    source_frag = group_connected_atoms(source_non_core, source_ct)
    dest_frag = group_connected_atoms(dest_non_core, dest_ct)

    # dictionary of attachment points, related fragments in the other molecule
    s_frag_visited = {}
    d_frag_visited = {}
    for (source_idx, dest_idx) in zip(source_core_atoms, dest_core_atoms):

        s_conn_atoms = find_connected_fragments(source_idx, source_frag, source_ct)
        d_conn_atoms = find_connected_fragments(dest_idx, dest_frag, dest_ct)

        # need to group close attachment points together
        s_frag_grouped = set()
        d_frag_grouped = set()
        d_frag_atom_grouped = set()

        for (s_frag_id, s_atoms) in s_conn_atoms:
            for s_atom_idx in s_atoms:
                dist_frag = []
                for (d_frag_id, d_atoms) in d_conn_atoms:
                    dist_sd = []
                    for d_atom_idx in d_atoms:
                        if not d_atom_idx in d_frag_atom_grouped:
                            dist_sd.append(
                                [
                                    numpy.linalg.norm(
                                        numpy.array(source_ct.Atoms[s_atom_idx].coor)
                                        - numpy.array(dest_ct.Atoms[d_atom_idx].coor)
                                    ),
                                    d_atom_idx,
                                ]
                            )
                    if dist_sd:
                        # find the d_fragment with mimimum distance and use the mimimun distance atom for the matched connection points
                        (dist, d_atom_idx) = sorted(dist_sd, key=lambda x: x[0])[0]
                        dist_frag.append((dist, d_frag_id, d_atoms, d_atom_idx))
                if dist_frag:
                    (dist, d_frag_id, d_atoms, d_atom_idx) = sorted(dist_frag, key=lambda x: x[0])[0]
                    d_frag_grouped.add(d_frag_id)
                    s_frag_grouped.add(s_frag_id)
                    d_frag_atom_grouped.add(d_atom_idx)
                    update_visted_fragments(s_frag_visited, s_frag_id, d_frag_id, source_idx, s_atoms)
                    update_visted_fragments(d_frag_visited, d_frag_id, s_frag_id, dest_idx, d_atoms)

        for (s_frag_id, s_atoms) in s_conn_atoms:
            if not s_frag_id in s_frag_grouped:
                update_visted_fragments(s_frag_visited, s_frag_id, -1, source_idx, s_atoms)

        for (d_frag_id, d_atoms) in d_conn_atoms:
            if not d_frag_id in d_frag_grouped:
                update_visted_fragments(d_frag_visited, d_frag_id, -1, dest_idx, d_atoms)

    frag_groups = []

    d_frag_grouped = {}
    s_frag_grouped = {}

    core_match = set([x for x in zip(source_core_atoms, dest_core_atoms)])

    source_att_point_in_frag_group = {}
    dest_att_point_in_frag_group = {}

    for s_frag in list(s_frag_visited):
        if s_frag in s_frag_grouped:
            continue

        (r_source_frags, r_dest_frags) = find_all_related_fragments(s_frag, s_frag_visited, d_frag_visited)

        f_group = FragmentGroup()
        f_group_id = len(frag_groups)

        for s_frag_idx in r_source_frags:
            s_frag_grouped[s_frag_idx] = f_group_id
            f_group.addSourceAtoms(source_frag[int(s_frag_idx)])
            s_att_points = s_frag_visited[s_frag_idx][1]
            d_frag_related_to_s = s_frag_visited[s_frag_idx][0]

            # no bridge atoms in dest for source fragment
            if not d_frag_related_to_s:
                for s_att in list(s_att_points):
                    att_points = generate_connection_points(int(s_att), s_att_points[s_att], None, [])
                    f_group.addAttachmentPoints(att_points)
                    register_att_point_fragmentgroup(
                        source_att_point_in_frag_group,
                        int(s_att),
                        s_att_points[s_att],
                        f_group_id,
                        f_group.getNumSourceFrags() - 1,
                    )
            # source fragment is turned to some other groups
            # try to pair attachment points first, then deal with surplus
            for d_frag in d_frag_related_to_s:
                f_group.addDestAtoms(dest_frag[int(d_frag)])
                d_frag_grouped[str(d_frag)] = f_group_id

                d_att_points = d_frag_visited[str(d_frag)][1]

                for s_att in list(s_att_points):
                    s_paired = False
                    for d_att in list(d_att_points):
                        if (int(s_att), int(d_att)) in core_match:
                            att_points = generate_connection_points(
                                int(s_att), s_att_points[s_att], int(d_att), d_att_points[d_att]
                            )
                            f_group.addAttachmentPoints(att_points)
                            register_att_point_fragmentgroup(
                                dest_att_point_in_frag_group,
                                int(d_att),
                                d_att_points[d_att],
                                f_group_id,
                                f_group.getNumDestFrags() - 1,
                            )
                            s_paired = True
                        else:
                            att_points = generate_connection_points(None, [], int(d_att), d_att_points[d_att])
                            f_group.addAttachmentPoints(att_points)
                            register_att_point_fragmentgroup(
                                dest_att_point_in_frag_group,
                                int(d_att),
                                d_att_points[d_att],
                                f_group_id,
                                f_group.getNumDestFrags() - 1,
                            )
                    if not s_paired:
                        att_points = generate_connection_points(int(s_att), s_att_points[s_att], None, [])
                        f_group.addAttachmentPoints(att_points)
                    register_att_point_fragmentgroup(
                        source_att_point_in_frag_group,
                        int(s_att),
                        s_att_points[s_att],
                        f_group_id,
                        f_group.getNumSourceFrags() - 1,
                    )

        frag_groups.append(f_group)

    # all the rest dest fragment, not sharing attachment points with any source fragments
    for d_frag in set(d_frag_visited) - set(d_frag_grouped):
        d_att_points = d_frag_visited[str(d_frag)][1]
        f_group = FragmentGroup()
        f_group.addDestAtoms(dest_frag[int(d_frag)])
        for d_att in list(d_att_points):
            att_points = generate_connection_points(None, [], int(d_att), d_att_points[d_att])
            f_group.addAttachmentPoints(att_points)
            register_att_point_fragmentgroup(
                dest_att_point_in_frag_group,
                int(d_att),
                d_att_points[d_att],
                len(frag_groups),
                f_group.getNumDestFrags() - 1,
            )

        frag_groups.append(f_group)

    for frag_idx in set(range(len(dest_frag))) - set(map(int, list(d_frag_visited))):
        f_group = FragmentGroup()
        f_group.addDestAtoms(dest_frag[frag_idx])
        f_group.addAttpointToDestFragment(0, (None, None))
        frag_groups.append(f_group)

    for frag_idx in set(range(len(source_frag))) - set(map(int, list(s_frag_visited))):
        f_group = FragmentGroup()
        f_group.addSourceAtoms(source_frag[frag_idx])
        f_group.addAttpointToSourceFragment(0, (None, None))
        frag_groups.append(f_group)

    for att_point in list(source_att_point_in_frag_group):
        (g_id, frag_id) = source_att_point_in_frag_group[att_point]
        frag_groups[g_id].addAttpointToSourceFragment(frag_id, att_point)
    for att_point in list(dest_att_point_in_frag_group):
        (g_id, frag_id) = dest_att_point_in_frag_group[att_point]
        frag_groups[g_id].addAttpointToDestFragment(frag_id, att_point)
    return (frag_groups, source_att_point_in_frag_group, dest_att_point_in_frag_group)
