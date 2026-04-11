import itertools

import networkx as nx

from ....utils import logger
from .fragment_group import get_mutating_fragment_groups


def get_structure_bonds_with_order(st):
    bond_res = []
    for bond in st.Bonds:
        atom = st.Atoms[bond.a1]
        for i in atom.connect:
            if bond.a2 == i:
                if atom.bond_type in ["s", "S", "bs", "es", "fs"]:
                    bond_res.append([bond.a1, bond.a2, 1])
                elif atom.bond_type in ["r", "D", "d", "bd", "ed", "fd"]:
                    bond_res.append([bond.a1, bond.a2, 2])
                elif atom.bond_type in ["t", "T", "ft"]:
                    bond_res.append([bond.a1, bond.a2, 3])
                else:
                    bond_res.append([bond.a1, bond.a2, 1])
    return bond_res


class FindDeletePath:
    def __init__(self, wt, mut, wt_core, mut_core):

        self.wt, self.mut = wt, mut
        self.wt_bonds = get_structure_bonds_with_order(wt)
        self.mut_bonds = get_structure_bonds_with_order(mut)

        self.wt_core = wt_core
        self.mut_core = mut_core
        self.fragment_groups, self.s_ab_in_g, self.d_ab_in_g = get_mutating_fragment_groups(wt_core, mut_core, wt, mut)
        self.wt_core_atom_groups = self.get_core_atom_groups(self.wt_bonds, self.wt_core)
        self.mut_core_atom_groups = self.get_core_atom_groups(self.mut_bonds, self.mut_core)

    def get_core_atom_groups(self, bonds, core_atoms):
        g = nx.Graph()
        g.add_nodes_from(core_atoms)
        core = set(core_atoms)
        for a, b, _ in bonds:
            if a in core and b in core:
                g.add_edge(a, b)
        return list(nx.connected_components(g))

    def find_delete_bond(self):
        wt_path = {}
        mut_path = {}
        for fg in self.fragment_groups:
            path_for_fg = self.find_bonds_to_delete(
                fg, self.wt_bonds, self.mut_bonds, self.wt_core_atom_groups, self.mut_core_atom_groups
            )
            wt_path.update(path_for_fg[0])
            mut_path.update(path_for_fg[1])
        self.wt_del, self.mut_del = wt_path, mut_path

    def find_bonds_to_delete(self, frag_group, wt_bonds, mut_bonds, wt_core_atom_groups, mut_core_atom_groups):
        wt_path = {}
        mut_path = {}

        def convert_fragment_to_graph(bonds, atoms):
            g = nx.Graph()
            g.add_nodes_from(atoms)
            g.graph["dont_cut"] = set()
            for a, b, order in bonds:
                if a in atoms and b in atoms:
                    g.add_edge(a, b, order=order)
                    if order == 3:
                        g.graph["dont_cut"].add(a)
                        g.graph["dont_cut"].add(b)
            return g

        for (s_id, s_atoms) in enumerate(frag_group.getSourceAtoms()):
            att_p = self.sort_attachment_points(frag_group.getSourceFragmentAttpoints(s_id), wt_core_atom_groups)
            s_graph = convert_fragment_to_graph(wt_bonds, s_atoms)
            wt_frag = self.select_bonds_to_delete(att_p, s_graph)
            wt_path.update(wt_frag)
        for (d_id, d_atoms) in enumerate(frag_group.getDestAtoms()):
            att_p = self.sort_attachment_points(frag_group.getDestFragmentAttpoints(d_id), mut_core_atom_groups)
            d_graph = convert_fragment_to_graph(mut_bonds, d_atoms)
            mut_frag = self.select_bonds_to_delete(att_p, d_graph)
            mut_path.update(mut_frag)
        return wt_path, mut_path

    def sort_attachment_points(self, attachment_points, connected_core_atoms):
        sort_list = []
        for (ank, dum) in attachment_points:
            for core_group in connected_core_atoms:
                if ank in core_group:
                    sort_list.append([(ank, dum), len(core_group)])
        return sorted(sort_list, key=lambda x: (x[1], x[0][0]))

    def select_bonds_to_delete(self, att_points, g):
        def bond_in_cycles(bond, cycles):
            (a, b) = bond
            for c in cycles:
                if (a in c) and (b in c):
                    return True
            return False

        num_to_break = max(0, len(att_points) - 1)
        # dictionary to store two paths for each broken bonds
        ret_path = {}
        if num_to_break == 0:
            return ret_path
        cycles = nx.cycle_basis(g)
        att_remained = set([att for (att, rank) in att_points])
        rank_dict = {a[0]: a[1] for a in att_points}
        att_pair_to_sort = [[(p[0][0], p[1][0]), p[0][1] + p[1][1]] for p in itertools.combinations((att_points), 2)]
        for (att1, att2) in [p[0] for p in sorted(att_pair_to_sort, key=lambda x: x[1])]:
            if (att1 not in att_remained) or (att2 not in att_remained) or (not (nx.has_path(g, att1[1], att2[1]))):
                continue
            path = nx.shortest_path(g, att1[1], att2[1])
            bonds = []
            for idx in range(len(path) - 1):
                edge = (path[idx], path[idx + 1])
                if (
                    not bond_in_cycles(
                        edge,
                        cycles
                        # exlude triple-bond atoms and double bond to keep original conformations
                    )
                    and edge[1] not in g.graph["dont_cut"]
                    and edge[0] not in g.graph["dont_cut"]
                    and g.get_edge_data(*edge)["order"] != 2
                ):
                    bonds.append([idx, edge])
            if bonds:
                idx, bond = self._find_bond_in_path(bonds)
                # pi starts form att1 to broken bond
                pi = [att1[0]]
                [pi.append(path[x]) for x in range(idx + 1)]
                # pj starts form att2 to broken bond
                pj = [path[x] for x in range(idx + 1, len(path))]
                pj.append(att2[0])
                pj.reverse()

                (bi, bj) = bond
                if bi > bj:
                    ret_path[(bj, bi)] = [pj, pi]
                else:
                    ret_path[(bi, bj)] = [pi, pj]
                g.remove_edges_from([bond])

                att_delete = self.disconnect_att(g, att_remained)
                if att_delete:
                    att_remained.remove(att_delete)
            else:
                # delete att with lower rank here
                (rank1, rank2) = (rank_dict[att1], rank_dict[att2])
                att_del = att2
                att_keep = att1
                if rank1 < rank2:
                    att_del, att_keep = att_keep, att_del

                att_remained.remove(att_del)
                (bi, bj) = att_del
                if bi in g.graph["dont_cut"] or bj in g.graph["dont_cut"]:
                    logger.warning("broken bond: (%d, %d) should not be cleaved " % (bi, bj))
                pi = [bi, bj]
                pj = [att_keep[0]]
                [pj.append(x) for x in path]

                if bi > bj:
                    ret_path[(bj, bi)] = [pj, pi]
                else:
                    ret_path[(bi, bj)] = [pi, pj]

            if len(ret_path) == num_to_break:
                return ret_path
        raise RuntimeError

    def disconnect_att(self, g, att):
        for frag in nx.connected_components(g):
            n_att_in_frag = 0
            last_att = None
            for anchor, bridge in att:
                if bridge in frag:
                    n_att_in_frag += 1
                    last_att = (anchor, bridge)
            if n_att_in_frag == 1:
                return last_att
        return None

    def _find_bond_in_path(self, bonds):
        sorted_bonds = sorted(bonds, key=lambda x: x[0], reverse=True)
        middle = len(bonds) // 2
        return sorted_bonds[middle]
