import math
import operator
from functools import reduce
from itertools import repeat
from multiprocessing import Manager

from ....utils import logger
from ...dual_topology import num_soft_bonds


class LightEdge:
    def __init__(self, edge):
        self.structs = edge.structs
        self.name = edge.name
        self.node_names = [edge.nodes[0].name, edge.nodes[1].name]
        self.edge_cores = [edge.nodes[0].core, edge.nodes[1].core]
        self.atom_mapping = None
        self.similarity = None
        self.similarity_detail = None
        self.is_core_hopping = edge.is_core_hopping

        if edge.atom_mapping:
            n0, n1 = list(edge.atom_mapping)
            self.atom_mapping = {n0: edge.atom_mapping[n0], n1: edge.atom_mapping[n1]}
        else:
            self.atom_mapping = None

class Rule(object):
    def __init__(self, **kwargs):
        self._args = kwargs

    def __call__(self, edge):
        raise NotImplementedError("__call__ not implemented in subclass of Rule")


def exp_delta(delta):
    BETA = 0.1
    return math.exp(-BETA * delta)


def calc_in_place_rmsd(molecule0, atom_list0, molecule1, atom_list1):
    l = len(atom_list0)
    sum = 0
    for i in range(l):
        distance = 0
        for j in [0, 1, 2]:
            distance += (molecule0.Atoms[atom_list0[i]].coor[j] - molecule1.Atoms[atom_list1[i]].coor[j]) ** 2
        sum += distance
    return pow(sum / l, 0.5)


class MCS(Rule):
    def delta(self, edge):
        struct0, struct1 = edge.structs
        atom_mapping_result = edge.atom_mapping
        n0_heavy = [atom.No for atom in struct0.Atoms if atom.atom_number > 1]
        n1_heavy = [atom.No for atom in struct1.Atoms if atom.atom_number > 1]
        atom_map_set = set(atom_mapping_result[edge.node_names[0]])
        n0_heavy_set = set(n0_heavy)
        intersection_set = n0_heavy_set & atom_map_set
        return len(n0_heavy) + len(n1_heavy) - 2 * len(intersection_set)

    def __call__(self, edge):
        return exp_delta(self.delta(edge))


class Charge(Rule):
    def __call__(self, edge):
        struct0, struct1 = edge.structs
        charge0 = sum(struct0.formal_charge)
        charge1 = sum(struct1.formal_charge)
        diff = int(math.fabs(charge0 - charge1))

        # charge hopping
        s0_map_atom, s1_map_atom = edge.atom_mapping.values()
        for i, formal_charge in enumerate(struct0.formal_charge):
            if formal_charge != 0 and i not in s0_map_atom:
                logger.debug(f"formal charge in dummy atoms: {struct0.name}, {edge.name}, {charge0}, {charge1}")
                diff = abs(formal_charge)
        for i, formal_charge in enumerate(struct1.formal_charge):
            if formal_charge != 0 and i not in s1_map_atom:
                logger.debug(f"formal charge in dummy atoms: {struct1.name}, {edge.name}, {charge0}, {charge1}")
                diff = max(diff, abs(formal_charge))
        return 0.1**diff


class MinimumNumberOfAtom(Rule):
    def __init__(self, cutoff=4):
        super().__init__(cutoff=cutoff)

    def __call__(self, edge):
        struct0, struct1 = edge.structs
        atom_mapping_result = edge.atom_mapping
        n0_heavy = [atom.No for atom in struct0.Atoms if atom.atom_number > 1]
        # n1_heavy = [atom.No for atom in struct1.Atoms if atom.atom_number > 1]
        atom_map_set = set(atom_mapping_result[edge.node_names[0]])
        n0_heavy_set = set(n0_heavy)
        intersection_set = n0_heavy_set & atom_map_set
        return int(len(intersection_set) > self._args["cutoff"])


class Rmsd(Rule):
    """
    Calculate the RSMD between the two molecules.

    This rule has two parameters:
    - switching  RMSD values greater than this will result in similarity score
                 close to zero.
    - steepness  This controls how fast the similarity score is switched between
                 1 and 0.
    """

    def __init__(self, switching=1.5, steepness=30.0):
        super().__init__(switching=switching, steepness=steepness)

    def __call__(self, edge):
        ct0, ct1 = edge.structs
        amap = edge.atom_mapping
        mcs0, mcs1 = amap[edge.node_names[0]], amap[edge.node_names[1]]
        if not len(mcs0):
            return 0
        else:
            rmsd = calc_in_place_rmsd(ct0, mcs0, ct1, mcs1)

        return (
            1
            if rmsd == 0
            else math.atan(self._args["steepness"] * (math.log(self._args["switching"]) - math.log(rmsd))) / math.pi
            + 0.5
        )


class Cutoff(Rule):
    """
    Check if the similarity score is greater or equal to a predefined cutoff
    value.
    Similarity score:
      = 0, if not
      = 1, if the similarity score is greater than equal to the cutoff value.
    """

    def __init__(self, cutoff=0.6):
        super().__init__(cutoff=cutoff)
        self._log = math.log(cutoff or 1e-10)

    def __call__(self, simi):
        return int(simi >= self._param["cutoff"])

    def similog(self, simi):
        return math.log(simi) / self._log


class SoftBond(Rule):
    def __call__(self, edge):
        ct0, ct1 = edge.structs
        amap = edge.atom_mapping
        mcs0, mcs1 = amap[edge.node_names[0]], amap[edge.node_names[1]]
        num_bonds = num_soft_bonds(ct0, ct1, mcs0, mcs1)
        return exp_delta(num_bonds)


def similarity(edge, result):
    score = {}
    score[MCS.__name__] = MCS()(edge)
    score[Charge.__name__] = Charge()(edge)
    score[MinimumNumberOfAtom.__name__] = MinimumNumberOfAtom()(edge)
    score[Rmsd.__name__] = Rmsd()(edge)
    score[SoftBond.__name__] = SoftBond()(edge)
    edge.similarity = round(reduce(operator.mul, score.values()), 4)
    edge.similarity_detail = score
    result[edge.name] = edge


def assign_similarity(edges, num_procs=1):
    #logger.info("Start calculating similarity ...")
    light_edges = [LightEdge(edge) for edge in edges]
    with Manager() as manager:
        manager_dict = manager.dict()
        with manager.Pool(processes=num_procs) as pool:
            pool.starmap(similarity, zip(light_edges, repeat(manager_dict, len(light_edges))))
            for edge in edges:
                edge.similarity = manager_dict[edge.name].similarity
                edge.similarity_detail = manager_dict[edge.name].similarity_detail
                edge.is_core_hopping = edge.similarity_detail[SoftBond.__name__] < 1
                edge.is_charge_hopping = edge.similarity_detail[Charge.__name__] < 1

    for edge in edges:
        #logger.info(f"{edge.name}: {edge.similarity:.3f}")
        if edge.is_charge_hopping:
            logger.warning(f"{edge.name}-core_hopping: {edge.is_core_hopping}")
        #logger.debug(f"{edge.name}-detail: {edge.similarity_detail}")
        
