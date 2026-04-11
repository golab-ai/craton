from copy import deepcopy
from dataclasses import dataclass, field
from typing import List

from .find_delete_bond import FindDeletePath
from .match import (
    AttachmentBondTopology,
    GetRestrain,
    TopologyMatch,
)


@dataclass
class RingBreakingData:
    bond_formation: List = field(default_factory=list)
    bond_delete: List = field(default_factory=list)
    dihedral_restrain: List = field(default_factory=list)
    improper_restrain: List = field(default_factory=list)

    def __eq__(self, other):
        return (
            set(self.bond_formation) == set(other.bond_formation)
            and set(self.bond_delete) == set(other.bond_delete)
            and set(self.dihedral_restrain) == set(other.dihedral_restrain)
            and set(self.improper_restrain) == set(other.improper_restrain)
        )


class RingBreaking:
    def __init__(self, mol1, mol2, matches):
        self.mol1 = mol1
        self.mol2 = mol2
        self.matches = matches

    def run(self):
        # wt, mut = self.edge.structs
        wt_copy, mut_copy = deepcopy(self.mol1), deepcopy(self.mol2)
        wt_core = list(self.matches.keys())
        mut_core = list(self.matches.values())

        # Find delete path
        self.del_path = FindDeletePath(wt_copy, mut_copy, wt_core, mut_core)
        self.del_path.find_delete_bond()

        # Matching the topology
        self.top_match = TopologyMatch(self.del_path)
        self.top_match.run_matching()

        # Updating attachment bond topology
        attach_match = AttachmentBondTopology(self.top_match)
        attach_match.run_matching()

        # Find restrain for attachment bonds and cleavage atoms
        self.restrain = GetRestrain(self.mol1, self.mol2, self.top_match, attach_match)
        self.restrain.get_restrain()

    def get_ring_brekaing_parameters(self):
        self.run()
        self.wt_result = self._generate_result(
            self.del_path.wt_del,
            self.top_match.bond_match.mut_delete_core,
            self.top_match.mut_core_map,
            self.restrain.wt_restrain,
            self.restrain.wt_restrain_cleavage,
        )
        self.mut_result = self._generate_result(
            self.del_path.mut_del,
            self.top_match.bond_match.wt_delete_core,
            self.top_match.wt_core_map,
            self.restrain.mut_restrain,
            self.restrain.mut_restrain_cleavage,
        )
        # self.edge.set_data("core_hopping_data", (self.wt_result, self.mut_result))
        # setattr(self.edge, "core_hopping_data", (self.wt_result, self.mut_result))

    def _generate_result(self, del_path, delete_core, core_map, restrain, restrain_cleavage):
        result = RingBreakingData()
        result.bond_delete += list(del_path)

        for item in delete_core:
            atom1, atom2 = core_map[item[0]], core_map[item[1]]
            if atom1 > atom2:
                atom1, atom2 = atom2, atom1
            result.bond_formation.append((atom1, atom2))

        if restrain is not None:
            for _, dihedral, improper in restrain.values():
                if dihedral is not None:
                    result.dihedral_restrain.append(dihedral)
                if improper is not None:
                    result.improper_restrain.append(improper)

        dihedrals, impropers = restrain_cleavage
        if dihedrals:
            result.dihedral_restrain += dihedrals
        if impropers:
            result.improper_restrain += impropers
        return result


class RingBreakingTopology:

    def __init__(self, fep_topology):
        self._mol1 = fep_topology.mol1
        self._mol2 = fep_topology.mol2

    def get_bond_list(self, mol1, mol2, dummy1, dummy2):
        m1_pair12, m1_pair13, m1_pair14, m1_pair1n = mol1.get_dummy_1234n_pairs()

