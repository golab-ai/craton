import math
from typing import Dict, List, Tuple, Union

import pandas as pd

from ...utils import logger
from ..mapping.graph import Graph


def get_population_from_pka_file(pka_file, ph=7.0):
    """
    Args:
        pka_file: pka csv file which contains ligand name, pka value, and pka type (acid or base)
        ph: ph of the solution
        for example
            l1, l2, pka, pka_type

    Returns:
        population
    """
    pka_df = pd.read_csv(pka_file, skiprows=1, names=["ligand1", "ligand2", "pka", "pka_type"], skipinitialspace=True)
    population_dict: Dict[Tuple[str, ...], List[float]] = {}
    for _, row in pka_df.iterrows():
        ligand1, ligand2 = row["ligand1"], row["ligand2"]
        pop = math.pow(10, float(row["pka"]) - ph)
        ratio = pop / (pop + 1)
        if row["pka_type"] == "acid":
            population_dict[(ligand1, ligand2)] = [ratio, 1 - ratio]
        elif row["pka_type"] == "base":
            population_dict[(ligand1, ligand2)] = [1 - ratio, ratio]
        else:
            raise RuntimeError("pka type must be base or acid")
    return population_dict


class PkaTautomerCorrection:
    def __init__(
        self,
        graph: Graph,
        pka_list: Union[List[Tuple[Tuple[str, ...], float]]] = None,
        population_list: Union[List[Tuple[Tuple[str, ...], List[float]]]] = None,
        pka_file: str = None,
        ph=7.0,
    ):
        self.graph = graph
        if pka_list is not None:
            self.population = self.get_population_from_pka(pka_list)
        elif population_list is not None:
            self.population = self.get_population_from_pka(pka_list)
        elif pka_file is not None:
            self.population = get_population_from_pka_file(pka_file, ph=ph)
        else:
            raise RuntimeError("must supply pka list or population list!")
        self.corrected_dg = self.get_corrected_dg()

    def get_population_from_pka(
        self, pka_list: List[Tuple[Tuple[str, ...], float]], ph=7.0
    ) -> Dict[Tuple[str, ...], List[float]]:
        node = self.graph.nodes_dict
        population_dict: Dict[Tuple[str, ...], List[float]] = {}
        for pair in pka_list:
            pop = math.pow(10, pair[1] - ph)
            ratio = pop / (pop + 1)
            left, right = pair[0][0], pair[0][1]
            if (left not in node) or (right not in node):
                logger.warning(f"cannot find {left} or {right} in graph, skip this pka")
                continue
            left_charge, right_charge = sum(node[left].struct.formal_charge), sum(node[right].struct.formal_charge)
            if left_charge - right_charge == 1:  # acid
                population_dict[pair[0]] = [ratio, 1 - ratio]
            else:  # base
                population_dict[pair[0]] = [1 - ratio, ratio]
        return population_dict

    def get_population_from_population_list(
        self, population_list: List[Tuple[Tuple[str, ...], List[float]]]
    ) -> Dict[Tuple[str, ...], List[float]]:
        population_dict: Dict[Tuple[str, ...], List[float]] = {}
        node = self.graph.nodes_dict
        for pair in population_list:
            for name in pair:
                if name not in node:
                    logger.warning(f"cannot find {name} in graph, skip!")
                    continue
            population_dict[pair[0]] = pair[1]
        return population_dict

    def get_corrected_dg(self):
        R = 1.9872041e-3
        T = 298.15
        RT = R * T
        populations_pairs = [(titles, populations) for titles, populations in self.population.items()]
        nodes = self.graph.nodes_dict
        corrected_dg = {lig: 0.0 for lig in nodes}
        for title, population in populations_pairs:
            for idx_i, title_i in enumerate(title):
                summation = 0.0
                pi = population[idx_i]
                for idx_j, title_j in enumerate(title):
                    if title_i != title_j:
                        i_dg = nodes[title_i].get_data("cc_dg")
                        j_dg = nodes[title_j].get_data("cc_dg")
                        if i_dg is None or j_dg is None:
                            logger.warning("cc dg not in graph, maybe cycle closure goes wrong, skip!")
                            continue
                        ddg = j_dg - i_dg
                        pj = population[idx_j]
                        summation += (pj / pi) * math.exp(-ddg / RT)
                pi_complex = 1 / (1 + summation)
                corrected_dg[title_i] = -RT * math.log(pi) + RT * math.log(pi_complex)
        return corrected_dg
