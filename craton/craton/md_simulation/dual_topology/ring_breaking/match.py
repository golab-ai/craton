from collections import defaultdict

from numpy import rad2deg

from .find_delete_bond import FindDeletePath


def linear_bonds(atoms):
    retval = set()
    for i in range(len(atoms) - 1):
        retval.add(frozenset((atoms[i], atoms[i + 1])))
    return retval


# for calculate similarity
def get_num_soft_bonds(wt, mut, wt_core, mut_core):
    del_path_obj = FindDeletePath(wt, mut, wt_core, mut_core)
    del_path_obj.find_delete_bond()
    wt_bonds = [(bond.a1, bond.a2) for bond in del_path_obj.wt.Bonds]
    mut_bonds = [(bond.a1, bond.a2) for bond in del_path_obj.mut.Bonds]
    top = TopologyMatch(del_path_obj=del_path_obj)
    pair_match = top.match_pairs(wt_bonds, mut_bonds)
    return (
        len(pair_match.wt_delete_core)
        + len(pair_match.mut_delete_core)
        + len(del_path_obj.wt_del)
        + len(del_path_obj.mut_del)
    )


class MatchedInteractions:
    def __init__(self, matched, wt_unmatched, wt_delete_core, wt_delete, mut_unmatched, mut_delete_core, mut_delete):
        # list of matched interaction pairs
        self.matched = frozenset([x for x in matched])
        # unmatched interactions excluding the ones in the deleted sets
        self.wt_unmatched = frozenset([x for x in wt_unmatched])
        self.mut_unmatched = frozenset([x for x in mut_unmatched])
        # core interaction deleted because of missing core bond
        self.wt_delete_core = frozenset([x for x in wt_delete_core])
        self.mut_delete_core = frozenset([x for x in mut_delete_core])
        # interaction involving dummpy atoms deleted
        # because of deleting bonds in non-interacting state
        self.wt_delete = frozenset([x for x in wt_delete])
        self.mut_delete = frozenset([x for x in mut_delete])


class TopologyMatch:
    def __init__(self, del_path_obj):
        self.del_path_obj = del_path_obj
        self.wt_core_map = {i: j for i, j in zip(del_path_obj.wt_core, del_path_obj.mut_core)}
        self.mut_core_map = {j: i for i, j in zip(del_path_obj.wt_core, del_path_obj.mut_core)}
        self.bond_match = None
        self.angle_match = None

    def find_terms_contain_bonds(self, terms, bond_in_term, bonds):
        deleted_terms = set()
        for b in bonds:
            deleted_terms.update(bond_in_term(terms, b))
        return deleted_terms

    def exclude_deleted_bonds_from_terms(self, terms, bond_in_term, dummy_bond_delete, core_bond_delete):
        term_dummy_del = self.find_terms_contain_bonds(terms, bond_in_term, dummy_bond_delete)
        term_core_del = self.find_terms_contain_bonds(terms - term_dummy_del, bond_in_term, core_bond_delete)
        return (terms - term_dummy_del - term_core_del, term_dummy_del, term_core_del)

    def match_pairs(self, wt_pairs, mut_pairs):
        matched_wt = []
        matched_mut = []
        wt_unmatched = set()
        mut_unmatched = set()
        wt_delete_core = set()
        mut_delete_core = set()
        for (wt_a, wt_b) in wt_pairs:
            if wt_a in self.wt_core_map:
                mut_a = self.wt_core_map[wt_a]
            else:
                wt_unmatched.add((wt_a, wt_b))
                continue
            if wt_b in self.wt_core_map:
                mut_b = self.wt_core_map[wt_b]
            else:
                wt_unmatched.add((wt_a, wt_b))
                continue

            if mut_b < mut_a:
                (mut_a, mut_b) = (mut_b, mut_a)

            if (mut_a, mut_b) in mut_pairs:
                matched_wt.append((wt_a, wt_b))
                matched_mut.append((mut_a, mut_b))
            else:
                wt_delete_core.add((wt_a, wt_b))

        for (mut_a, mut_b) in set(mut_pairs) - set(matched_mut):
            if (mut_a in self.mut_core_map) and (mut_b in self.mut_core_map):
                mut_delete_core.add((mut_a, mut_b))
            else:
                mut_unmatched.add((mut_a, mut_b))
        return MatchedInteractions(
            list(zip(matched_wt, matched_mut)),
            wt_unmatched,
            wt_delete_core,
            set([]),
            mut_unmatched,
            mut_delete_core,
            set([]),
        )

    def match_14pairs(self):
        self.del_path_obj.wt.create_intra_nonbond()
        self.del_path_obj.mut.create_intra_nonbond()
        wt_14pair = [(pair.a1, pair.a2) for pair in self.del_path_obj.wt.Pair14]
        mut_14pair = [(pair.a1, pair.a2) for pair in self.del_path_obj.mut.Pair14]
        self.pair14_match = self.match_pairs(wt_14pair, mut_14pair)

    def match_exclusion(self):
        wt_12pair = set([(pair.a1, pair.a2) for pair in self.del_path_obj.wt.Pair12])
        mut_12pair = set([(pair.a1, pair.a2) for pair in self.del_path_obj.mut.Pair12])
        wt_13pair = set([(pair.a1, pair.a2) for pair in self.del_path_obj.wt.Pair13])
        mut_13pair = set([(pair.a1, pair.a2) for pair in self.del_path_obj.mut.Pair13])
        wt_14pair = set([(pair.a1, pair.a2) for pair in self.del_path_obj.wt.Pair14])
        mut_14pair = set([(pair.a1, pair.a2) for pair in self.del_path_obj.mut.Pair14])
        wt_excl = wt_12pair | wt_13pair | wt_14pair
        mut_excl = mut_12pair | mut_13pair | mut_14pair
        self.exclusion_match = self.match_pairs(wt_excl, mut_excl)

    def match_bonds(self):
        wt_bonds = [(bond.a1, bond.a2) for bond in self.del_path_obj.wt.Bonds]
        mut_bonds = [(bond.a1, bond.a2) for bond in self.del_path_obj.mut.Bonds]
        wt_bonds_delete = list(self.del_path_obj.wt_del)
        mut_bonds_delete = list(self.del_path_obj.mut_del)

        pm = self.match_pairs(wt_bonds, mut_bonds)
        wt_delete = pm.wt_delete.union(wt_bonds_delete)
        mut_delete = pm.mut_delete.union(mut_bonds_delete)
        self.bond_match = MatchedInteractions(
            pm.matched,
            pm.wt_unmatched - wt_delete,
            pm.wt_delete_core,
            wt_delete,
            pm.mut_unmatched - mut_delete,
            pm.mut_delete_core,
            mut_delete,
        )

    def match_angles(self):
        def find_angles_contain_bond(angles, bond):
            ret_val = set()
            for angle in angles:
                if frozenset(bond) in linear_bonds(angle):
                    ret_val.add(angle)
            return ret_val

        wt_angles = [(angle.a1, angle.a2, angle.a3) for angle in self.del_path_obj.wt.Angles]
        mut_angles = [(angle.a1, angle.a2, angle.a3) for angle in self.del_path_obj.mut.Angles]

        matched_wt = []
        matched_mut = []
        wt_unmatched = set()
        mut_unmatched = set()
        wt_delete = set()
        mut_delete = set()
        wt_delete_core = set()
        mut_delete_core = set()
        for (wt_i, wt_j, wt_k) in wt_angles:
            if wt_i in self.wt_core_map:
                mut_i = self.wt_core_map[wt_i]
            else:
                wt_unmatched.add((wt_i, wt_j, wt_k))
                continue
            if wt_j in self.wt_core_map:
                mut_j = self.wt_core_map[wt_j]
            else:
                wt_unmatched.add((wt_i, wt_j, wt_k))
                continue
            if wt_k in self.wt_core_map:
                mut_k = self.wt_core_map[wt_k]
            else:
                wt_unmatched.add((wt_i, wt_j, wt_k))
                continue

            if mut_i > mut_k:
                (mut_i, mut_k) = (mut_k, mut_i)

            if (mut_i, mut_j, mut_k) in mut_angles:
                matched_wt.append((wt_i, wt_j, wt_k))
                matched_mut.append((mut_i, mut_j, mut_k))
            else:
                wt_delete_core.add((wt_i, wt_j, wt_k))

        for (mut_i, mut_j, mut_k) in set(mut_angles) - set(matched_mut):
            if mut_i in self.mut_core_map and mut_k in self.mut_core_map and mut_j in self.mut_core_map:
                mut_delete_core.add((mut_i, mut_j, mut_k))
            else:
                mut_unmatched.add((mut_i, mut_j, mut_k))

        (wt_final_unmatched, wt_delete, wt_more_delete_core) = self.exclude_deleted_bonds_from_terms(
            wt_unmatched, find_angles_contain_bond, self.bond_match.wt_delete, self.bond_match.wt_delete_core
        )
        wt_delete_core.update(wt_more_delete_core)
        (mut_final_unmatched, mut_delete, mut_more_delete_core) = self.exclude_deleted_bonds_from_terms(
            mut_unmatched, find_angles_contain_bond, self.bond_match.mut_delete, self.bond_match.mut_delete_core
        )
        mut_delete_core.update(mut_more_delete_core)

        self.angle_match = MatchedInteractions(
            list(zip(matched_wt, matched_mut)),
            wt_final_unmatched,
            wt_delete_core,
            wt_delete,
            mut_final_unmatched,
            mut_delete_core,
            mut_delete,
        )

    def match_dihedrals(self):
        def find_dihedrals_contain_bond(dihedrals, bond):
            ret_val = set()
            for dihedral in dihedrals:
                if frozenset(bond) in linear_bonds(dihedral):
                    ret_val.add(dihedral)
            return ret_val

        wt_dihedrals = [
            (dihedral.a1, dihedral.a2, dihedral.a3, dihedral.a4) for dihedral in self.del_path_obj.wt.Dihedrals
        ]
        mut_dihedrals = [
            (dihedral.a1, dihedral.a2, dihedral.a3, dihedral.a4) for dihedral in self.del_path_obj.mut.Dihedrals
        ]
        matched_wt = []
        matched_mut = []
        wt_unmatched = set()
        mut_unmatched = set()
        wt_delete = set()
        mut_delete = set()
        wt_delete_core = set()
        mut_delete_core = set()

        for (wt_i, wt_j, wt_k, wt_l) in wt_dihedrals:
            if wt_i in self.wt_core_map:
                mut_i = self.wt_core_map[wt_i]
            else:
                wt_unmatched.add((wt_i, wt_j, wt_k, wt_l))
                continue
            if wt_j in self.wt_core_map:
                mut_j = self.wt_core_map[wt_j]
            else:
                wt_unmatched.add((wt_i, wt_j, wt_k, wt_l))
                continue
            if wt_k in self.wt_core_map:
                mut_k = self.wt_core_map[wt_k]
            else:
                wt_unmatched.add((wt_i, wt_j, wt_k, wt_l))
                continue
            if wt_l in self.wt_core_map:
                mut_l = self.wt_core_map[wt_l]
            else:
                wt_unmatched.add((wt_i, wt_j, wt_k, wt_l))
                continue

            if mut_i > mut_l:  # the diherdal index always have i < l
                (mut_i, mut_j, mut_k, mut_l) = (mut_l, mut_k, mut_j, mut_i)

            if (mut_i, mut_j, mut_k, mut_l) in mut_dihedrals:
                matched_wt.append((wt_i, wt_j, wt_k, wt_l))
                matched_mut.append((mut_i, mut_j, mut_k, mut_l))
            else:
                wt_delete_core.add((wt_i, wt_j, wt_k, wt_l))

        for (mut_i, mut_j, mut_k, mut_l) in set(mut_dihedrals) - set(matched_mut):
            if (
                mut_i in self.mut_core_map
                and mut_k in self.mut_core_map
                and mut_j in self.mut_core_map
                and mut_l in self.mut_core_map
            ):
                mut_delete_core.add((mut_i, mut_j, mut_k, mut_l))
            else:
                mut_unmatched.add((mut_i, mut_j, mut_k, mut_l))

        (wt_final_unmatched, wt_delete, wt_more_delete_core) = self.exclude_deleted_bonds_from_terms(
            wt_unmatched, find_dihedrals_contain_bond, self.bond_match.wt_delete, self.bond_match.wt_delete_core
        )
        wt_delete_core.update(wt_more_delete_core)
        (mut_final_unmatched, mut_delete, mut_more_delete_core) = self.exclude_deleted_bonds_from_terms(
            mut_unmatched, find_dihedrals_contain_bond, self.bond_match.mut_delete, self.bond_match.mut_delete_core
        )
        mut_delete_core.update(mut_more_delete_core)

        self.dihedral_match = MatchedInteractions(
            list(zip(matched_wt, matched_mut)),
            wt_final_unmatched,
            wt_delete_core,
            wt_delete,
            mut_final_unmatched,
            mut_delete_core,
            mut_delete,
        )

    def match_impropers(self):
        def bonded(i, j, k, l, bonds):
            return (
                ((i, j) in bonds or (j, i) in bonds)
                and ((i, k) in bonds or (k, i) in bonds)
                and ((i, l) in bonds or (l, i) in bonds)
            )

        def find_impropers_contain_bond(impropers, bond):
            ret_val = set()
            for improper in impropers:
                (i, j, k, l) = improper
                if frozenset(bond) in set(
                    (
                        frozenset((i, k)),
                        frozenset((i, j)),
                        frozenset((i, l)),
                        frozenset((k, i)),
                        frozenset((l, i)),
                        frozenset((j, i)),
                    )
                ):
                    ret_val.add(improper)
            return ret_val
        if not hasattr(self.del_path_obj.wt, "Impropers") or not hasattr(self.del_path_obj.mut, "Impropers"):
            self.improper_match = None
            return

        wt_impropers = set(
            [(improper.a1, improper.a2, improper.a3, improper.a4) for improper in self.del_path_obj.wt.Impropers]
        )
        mut_impropers = set(
            [(improper.a1, improper.a2, improper.a3, improper.a4) for improper in self.del_path_obj.mut.Impropers]
        )
        mut_bonds = set([(bond.a1, bond.a2) for bond in self.del_path_obj.wt.Bonds])

        matched_wt = []
        matched_mut = []
        wt_unmatched = set()
        mut_unmatched = set()
        wt_delete = set()
        mut_delete = set()
        wt_delete_core = set()
        mut_delete_core = set()

        for (wt_i, wt_j, wt_k, wt_l) in wt_impropers:
            if wt_i in self.wt_core_map:
                mut_i = self.wt_core_map[wt_i]
            else:
                wt_unmatched.add((wt_i, wt_j, wt_k, wt_l))
                continue
            if wt_j in self.wt_core_map:
                mut_j = self.wt_core_map[wt_j]
            else:
                wt_unmatched.add((wt_i, wt_j, wt_k, wt_l))
                continue
            if wt_k in self.wt_core_map:
                mut_k = self.wt_core_map[wt_k]
            else:
                wt_unmatched.add((wt_i, wt_j, wt_k, wt_l))
                continue
            if wt_l in self.wt_core_map:
                mut_l = self.wt_core_map[wt_l]
            else:
                wt_unmatched.add((wt_i, wt_j, wt_k, wt_l))
                continue

            # the i is the center atom
            if bonded(mut_i, mut_j, mut_k, mut_l, mut_bonds):
                if (mut_i, mut_j, mut_k, mut_l) in mut_impropers:
                    matched_wt.append((wt_i, wt_j, wt_l, wt_l))
                    matched_mut.append((mut_i, mut_j, mut_k, mut_l))
                if (mut_i, mut_k, mut_j, mut_l) in mut_impropers:
                    matched_wt.append((wt_i, wt_j, wt_l, wt_l))
                    matched_mut.append((mut_i, mut_k, mut_j, mut_l))
            else:
                wt_delete_core.add((wt_i, wt_j, wt_k, wt_l))

        for (mut_i, mut_j, mut_k, mut_l) in set(mut_impropers) - set(matched_mut):
            if (
                mut_i in self.mut_core_map
                and mut_k in self.mut_core_map
                and mut_j in self.mut_core_map
                and mut_l in self.mut_core_map
            ):
                mut_delete_core.add((mut_i, mut_j, mut_k, mut_l))
            else:
                mut_unmatched.add((mut_i, mut_j, mut_k, mut_l))

        (wt_final_unmatched, wt_delete, wt_more_delete_core) = self.exclude_deleted_bonds_from_terms(
            wt_unmatched, find_impropers_contain_bond, self.bond_match.wt_delete, self.bond_match.wt_delete_core
        )
        wt_delete_core.update(wt_more_delete_core)
        (mut_final_unmatched, mut_delete, mut_more_delete_core) = self.exclude_deleted_bonds_from_terms(
            mut_unmatched, find_impropers_contain_bond, self.bond_match.mut_delete, self.bond_match.mut_delete_core
        )
        mut_delete_core.update(mut_more_delete_core)

        self.improper_match = MatchedInteractions(
            list(zip(matched_wt, matched_mut)),
            wt_final_unmatched,
            wt_delete_core,
            wt_delete,
            mut_final_unmatched,
            mut_delete_core,
            mut_delete,
        )

    def run_matching(self):
        self.match_bonds()
        self.match_14pairs()
        self.match_exclusion()
        self.match_angles()
        self.match_dihedrals()
        self.match_impropers()


class AttachmentBondTopology:
    def __init__(self, top_match):
        self.top_match: TopologyMatch = top_match
        self.wt_att_bonds = self.remove_deleted_att(
            self.top_match.del_path_obj.s_ab_in_g, self.top_match.bond_match.wt_delete
        )
        self.mut_att_bonds = self.remove_deleted_att(
            self.top_match.del_path_obj.d_ab_in_g, self.top_match.bond_match.mut_delete
        )
        self.att_dict = defaultdict(dict)
        self.angle_kept = {}
        self.dihedral_kept = {}
        self.improper_kept = {}

    def remove_deleted_att(self, att_dict, bond_delete):
        ret_val = set()
        for (i, j) in att_dict:
            if i > j:
                (i, j) = (j, i)
            if not (i, j) in bond_delete:
                ret_val.add((i, j))
        return ret_val

    def find_angle_for_att(self):
        def find_angles_contain_bond(angle_set, bond):
            ret_val = set()
            for angle in angle_set:
                if frozenset(bond) in linear_bonds(angle):
                    ret_val.add(angle)
            return ret_val

        for item in ["wt", "mut"]:
            if item == "wt":
                (unmatched, attach_bond, core) = (
                    self.top_match.angle_match.wt_unmatched,
                    self.wt_att_bonds,
                    list(self.top_match.wt_core_map.keys()),
                )
            elif item == "mut":
                (unmatched, attach_bond, core) = (
                    self.top_match.angle_match.mut_unmatched,
                    self.mut_att_bonds,
                    list(self.top_match.wt_core_map.values()),
                )
            angle_set = set([angle for angle in unmatched])
            for bond in attach_bond:
                core_angle = set()
                for (i, j, k) in find_angles_contain_bond(angle_set, bond):
                    if set((i, j)) == set(bond):
                        if k in core:
                            core_angle.add((i, j, k))
                    elif set((j, k)) == set(bond):
                        if i in core:
                            core_angle.add((i, j, k))
                angle_set -= core_angle  # two dummy, one core
                self.att_dict[item][bond] = [core_angle]
            self.angle_kept[item] = angle_set

    def find_dihedral_for_att(self):
        def find_dihedrals_contain_bond(dihedrals, bond):
            ret_val = set()
            for dihedral in dihedrals:
                if frozenset(bond) in linear_bonds(dihedral):
                    ret_val.add(dihedral)
            return ret_val

        for item in ["wt", "mut"]:
            if item == "wt":
                (unmatched, core) = (
                    self.top_match.dihedral_match.wt_unmatched,
                    list(self.top_match.wt_core_map.keys()),
                )
            elif item == "mut":
                (unmatched, core) = (
                    self.top_match.dihedral_match.mut_unmatched,
                    list(self.top_match.wt_core_map.values()),
                )

            um = set([x for x in unmatched])
            for bond in list(self.att_dict[item]):
                core_dihe = set()
                for (i, j, k, l) in find_dihedrals_contain_bond(um, bond):
                    if set((i, j)) == set(bond):
                        if k in core or l in core:
                            core_dihe.add((i, j, k, l))
                    elif set((j, k)) == set(bond):
                        if i in core or l in core:
                            core_dihe.add((i, j, k, l))
                    elif set((k, l)) == set(bond):
                        if i in core or j in core:
                            core_dihe.add((i, j, k, l))
                self.att_dict[item][bond].append(core_dihe)
                um -= core_dihe  # three dummy one core
            self.dihedral_kept[item] = um

    def find_impropers_for_att(self):

        if not self.top_match.improper_match:
            return

        def find_impropers_contain_bond(impropers, bond):
            ret_val = set()
            for improper in impropers:
                (i, j, k, l) = improper
                if frozenset(bond) in set((frozenset((i, k)), frozenset((i, j)), frozenset((i, l)))):
                    ret_val.add(improper)
            return ret_val

        for item in ["wt", "mut"]:
            if item == "wt":
                (unmatched, core) = (
                    self.top_match.improper_match.wt_unmatched,
                    list(self.top_match.wt_core_map.keys()),
                )
            elif item == "mut":
                (unmatched, core) = (
                    self.top_match.improper_match.mut_unmatched,
                    list(self.top_match.wt_core_map.values()),
                )
            um = set([x for x in unmatched])
            for bond in list(self.att_dict[item]):
                core_impt = set()
                for (i, j, k, l) in find_impropers_contain_bond(um, bond):
                    if set((i, j)) == set(bond):
                        if k in core or l in core:
                            core_impt.add((i, j, k, l))
                    elif set((i, k)) == set(bond):
                        if i in core or l in core:
                            core_impt.add((i, j, k, l))
                    elif set((i, l)) == set(bond):
                        if i in core or j in core:
                            core_impt.add((i, j, k, l))
                self.att_dict[item][bond].append(core_impt)
                um -= core_impt  # three dummy one core
            self.improper_kept[item] = um

    def run_matching(self):
        self.find_angle_for_att()
        self.find_dihedral_for_att()
        self.find_impropers_for_att()


class GetRestrain:
    def __init__(self, wt, mut, top_match, attach_match):
        self.wt = wt
        self.mut = mut
        self.wt_angles_param = {(angle.a1, angle.a2, angle.a3): angle.parameter for angle in wt.Angles}
        self.mut_angles_param = {(angle.a1, angle.a2, angle.a3): angle.parameter for angle in mut.Angles}
        self.top_match: TopologyMatch = top_match
        self.attach_match: AttachmentBondTopology = attach_match

    def match_angle(self, core_map, angle):
        (i, j, k) = angle
        if i in core_map and j in core_map and k in core_map:
            mut_i = core_map[i]
            mut_j = core_map[j]
            mut_k = core_map[k]
            if mut_i > mut_k:
                mut_i, mut_k = mut_k, mut_i
            return (mut_i, mut_j, mut_k)

        return (None, None, None)

    def get_restrain(self):
        self.wt_restrain = self._restrain_for_attachment_bonds(
            self.wt,
            self.wt_angles_param,
            self.mut_angles_param,
            self.attach_match.att_dict["wt"],
            self.top_match.wt_core_map,
            self.top_match.bond_match.wt_delete_core,
        )
        self.mut_restrain = self._restrain_for_attachment_bonds(
            self.mut,
            self.mut_angles_param,
            self.wt_angles_param,
            self.attach_match.att_dict["mut"],
            self.top_match.mut_core_map,
            self.top_match.bond_match.mut_delete_core,
        )
        self.wt_restrain_cleavage = self._restrain_for_dummy_at_bond_cleavage(
            self.wt, self.top_match.del_path_obj.wt_del, self.top_match.del_path_obj.wt_core
        )

        self.mut_restrain_cleavage = self._restrain_for_dummy_at_bond_cleavage(
            self.mut, self.top_match.del_path_obj.mut_del, self.top_match.del_path_obj.mut_core
        )

    def _restrain_for_attachment_bonds(self, wt, angle_param, other_angle_param, att_dict, core_map, delete_core_bonds):
        """
        Find all the possible atoms for dihedral or improper restraint across attachment bonds
        ct: mol object
        att_dict: dictionary of interactions keyed by attachment bonds
        core_map: map core atoms to the ones in the other molecule
        angles_param: parameter dictionary for angles in current molecule
        other_angles_param: parameter dictionary of angles for the other molecule involving
                            alchemical change
        deleted_core_bonds: frozenset of core-core bonds to be deleted
        return: dictionary: key = attachment-bond, value = a list of tuples, where each tuple
                is (angle, proper/None, improper/None). where the proper or improper is the
                is the candidate dihedral interation to be restrained.
        """

        def bond_not_in_set(bond, bond_set):
            return not (frozenset(bond) in [frozenset(b) for b in bond_set])

        ret_val = {}
        bond_visited = set()
        # group attachment points together
        att_grouped = {}
        for bond in att_dict:
            (ank, dum) = bond
            if dum in core_map:
                (dum, ank) = bond
            att_grouped.setdefault(ank, []).append(bond)
        for k in att_grouped:
            if len(att_grouped[k]) > 1:
                extra_core = None
                for bond in att_grouped[k]:
                    # sort the angles so that more stable atom gets picked first
                    angles_w_param = [(a, angle_param[a][1]) for a in att_dict[bond][0]]
                    angles_sorted = sorted(angles_w_param, key=lambda x: x[1], reverse=True)
                    for (angle, _) in angles_sorted:
                        more_core = set(angle) - set(bond)
                        if more_core:
                            extra_core = more_core.pop()
                            break
                if extra_core:
                    # pick the first dummy atom for a set of improper restraints
                    dum0 = (set(att_grouped[k][0]) - set([k])).pop()
                    bond_visited.add(att_grouped[k][0])
                    # no need to restrain the first dummy, it will be included later
                    angle = (dum0, k, extra_core) if dum0 < extra_core else (extra_core, k, dum0)
                    ret_val[att_grouped[k][0]] = (angle, None, None)

                    already_restrained = set([extra_core, dum0])

                    for bond in att_grouped[k][1:]:
                        dum = (set(bond) - set([k])).pop()
                        bond_visited.add(bond)
                        angle = (dum, k, extra_core) if dum < extra_core else (extra_core, k, dum)
                        ret_val[bond] = (angle, None, (k, dum0, extra_core, dum))
                        already_restrained.add(dum)
                    # leftover core atom, need one more restraint
                    core_atom_left = (
                        set(
                            [
                                a
                                for a in wt.Atoms[k].connect
                                if (a in core_map and bond_not_in_set((k, a), delete_core_bonds))
                            ]
                        )
                        - already_restrained
                    )
                    if core_atom_left:
                        c_atom = core_atom_left.pop()
                        (angle, dihe, impt) = ret_val[att_grouped[k][0]]
                        core_angle = (extra_core, k, c_atom) if extra_core < c_atom else (c_atom, k, extra_core)
                        matched_angle = self.match_angle(core_map, core_angle)

                        # need to check this angle and its match to avoid colinear
                        # improper
                        if (matched_angle in other_angle_param) and (core_angle in angle_param):
                            (a0, k0) = angle_param[core_angle]
                            (a1, k1) = other_angle_param[matched_angle]
                            if (abs(rad2deg(a0) - 180.0) > 0.8) and (abs(rad2deg(a1) - 180.0) > 0.8):
                                ret_val[att_grouped[k][0]] = (angle, None, (dum0, extra_core, k, c_atom))
        for bond in att_dict:
            if bond in bond_visited:
                continue
            try:
                (angles, dihes, impts) = att_dict[bond]
            except ValueError:
                (angles, dihes) = att_dict[bond]
            if len(angles) == 1:
                # special case for partially matched triple bond
                ret_val[bond] = (list(angles)[0], None, None)
                continue

            tmp_for_sort = []
            for angle in angles:
                (ka, a0) = angle_param[angle]
                if ka:
                    tmp_for_sort.append([ka, a0, angle])
            if len(tmp_for_sort):
                a_sort = sorted(tmp_for_sort, key=lambda x: x[1], reverse=True)
                (ka, a0, angle) = a_sort[0]
                (a_i, a_j, a_k) = angle

                # Avoids using the attachment bond as the line of intersection.
                # If the attachment bond is the (a_i, a_j), we choose improper dihedral to be
                # (a_i, a_j, a_k, l); if not, (a_k, a_j, a_i, l), where l is the fouth atom to be
                # found below.
                (impt_i, impt_j, impt_k) = (a_i, a_j, a_k) if set((a_i, a_j)) == set(bond) else (a_k, a_j, a_i)

                # Now we are to find the 'l' atom. We loop through all other angles that involve
                # the attachment bond.
                triplet = (impt_i, impt_j, impt_k)
                for (a0, ka, angle) in a_sort[1:]:
                    if abs(rad2deg(a0) - 180.0) < 0.5:
                        # We don't consider the angle if the three atoms are colinear.
                        continue

                    # Resets to the original triplet, because impt_* might be altered during the
                    # previous iteration.
                    impt_i, impt_j, impt_k = triplet
                    (a_i1, a_j1, a_k1) = angle

                    if set((a_j1, a_k1)) == set(bond):
                        impt_l = a_i1
                    else:
                        impt_l = a_k1

                    # Now we have an 'l' atom candidate, and it should be a core atom.
                    # To be a valid l atom, the j, k, and l atoms must not be colinear in the B state.
                    o_a = self.match_angle(core_map, (impt_l, impt_j, impt_k))

                    if o_a in other_angle_param:
                        (a0, ka) = other_angle_param[o_a]
                        if abs(rad2deg(a0) - 180.0) < 0.5:
                            # linear molecule, can not define improper angle, skipping
                            continue
                        else:
                            # found improper
                            if impt_i > impt_l:
                                (impt_i, impt_l) = (impt_l, impt_i)
                            ret_val[bond] = (
                                (a_i, a_j, a_k),
                                None,
                                (impt_j, impt_i, impt_k, impt_l),
                            )  # first atom is center atom
                            # Are we good with one candidate?
                            break

                if bond not in ret_val:
                    # need to look for a dihedral angle to restrain
                    for (a0, ka, angle) in a_sort:
                        (a_i, a_j, a_k) = angle
                        for (dihe_a1, dihe_a2, dihe_a3, dihe_a4) in dihes:
                            core_atom = (
                                (dihe_a1 in core_map)
                                + (dihe_a2 in core_map)
                                + (dihe_a3 in core_map)
                                + (dihe_a4 in core_map)
                            )
                            if core_atom == 3:
                                if set((a_i, a_j, a_k)) == set((dihe_a1, dihe_a2, dihe_a3)) or set(
                                    (a_i, a_j, a_k)
                                ) == set(dihe_a2, dihe_a3, dihe_a4):
                                    ret_val[bond] = ((a_i, a_j, a_k), (dihe_a1, dihe_a2, dihe_a3, dihe_a4), None)
                                    break
        return ret_val

    def _restrain_for_dummy_at_bond_cleavage(self, wt, path_dict, core):
        def path_contains_cleaved_bond(path, bonds):
            bonds_in_path = set()
            for i in range(len(path) - 1):
                bonds_in_path.add(frozenset([path[i], path[i + 1]]))
            for b in bonds:
                if b in bonds_in_path:
                    return True
            return False

        dummy_atoms = set(range(len(wt.Atoms))) - set(core)
        dihe = []
        impt = []

        all_path = set()
        cleaved_bonds = set()
        for path in path_dict:
            cleaved_bonds.add(frozenset(path))
            for p in path_dict[path]:
                all_path.update(p)
        visited = set()
        for k in path_dict:
            for atom, path in zip(k, path_dict[k]):
                if path_contains_cleaved_bond(path, cleaved_bonds):
                    continue
                terminal_dummy_atoms = [a for a in set(wt.Atoms[atom].connect) & set(dummy_atoms) - all_path]
                core_bonded_to_path = [a for a in set(wt.Atoms[path[0]].connect) & set(core)]
                (d0, imt0) = self._get_restraint_for_atom(
                    path, atom, core_bonded_to_path, terminal_dummy_atoms, visited
                )
                dihe.extend(d0)
                impt.extend(imt0)
        return (dihe, impt)

    def _get_restraint_for_atom(self, path, atom, core_bonded_to_path, terminal_dummy_atoms, visited):
        def get_restraint_from_path(path, dummy_atoms, core_atoms=()):

            natom_need = max(4 - len(dummy_atoms), 2)
            more_atoms = []
            if len(path) < natom_need:
                # path is not long enough, need to find one more core atom
                if core_atoms:
                    more_atoms.append(core_atoms[0])
                else:
                    # no more core atoms, cannot define dihedral angle
                    return []
                more_atoms.extend(path)
            else:
                more_atoms.extend(path[len(path) - natom_need :])

            n_dummy = len(dummy_atoms)
            if n_dummy == 1:
                if len(more_atoms) == 3:
                    return [(more_atoms[0], more_atoms[1], more_atoms[2], dummy_atoms[0])]
                else:
                    # cannot find enough core atoms for restraint
                    return []
            elif n_dummy >= 2:
                ret_val = []
                for i in range(len(dummy_atoms) - 1):
                    ret_val.append((dummy_atoms[i], more_atoms[0], more_atoms[1], dummy_atoms[i + 1]))
                return ret_val
            else:
                # no dummy atoms to restrain
                return []

        dihe = []
        impt = []

        if not atom in visited:
            visited.add(atom)
            if len(terminal_dummy_atoms) == 1:
                dihe.extend(get_restraint_from_path(path, terminal_dummy_atoms, core_bonded_to_path))
            elif len(terminal_dummy_atoms) > 1:
                impt.extend(get_restraint_from_path(path, terminal_dummy_atoms))
        return (dihe, impt)


def match_interactions(wt, mut, wt_core, mut_core):
    delete_path_obj = FindDeletePath(wt, mut, wt_core, mut_core)

    # Find delete path
    delete_path_obj.find_delete_bond()

    # Mathcing the topoloty
    top_match = TopologyMatch(delete_path_obj)
    top_match.run_matching()

    # Updating attachment bond topology
    attach_match = AttachmentBondTopology(top_match)
    attach_match.run_matching()

    # Find retrain for attachment bonds and cleavage atoms
    restrain = GetRestrain(wt, mut, top_match, attach_match)
    restrain.get_restrain()
