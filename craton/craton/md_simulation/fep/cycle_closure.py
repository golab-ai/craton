import heapq
import json
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Tuple

import networkx as nx
import numpy as np
import pandas as pd

from ...utils import logger
from ..mapping.graph import Graph, Node
from .pka_correction import PkaTautomerCorrection

TOLERANCE = 1e-9
EXPERIMENT_ERROR = 0.16
NODE_CUTOFF = 6
NODE_CUTOFF_STAR_GRAPH = 3

# for UI interface


def run_cycle_closure_from_json_result(result_str: str) -> str:
    result_json = json.loads(result_str)
    node_start = result_json["summary"]["pair"]["begin"]
    node_end = result_json["summary"]["pair"]["end"]
    ddg = result_json["summary"]["pair"]["ddG"]
    error = result_json["summary"]["pair"]["error"]
    if result_json["summary"].get("exp") is None or result_json["summary"]["exp"].get("name") is None:
        exp_data = None
    else:
        exp_data = {
            key: value for key, value in zip(result_json["summary"]["exp"]["name"], result_json["summary"]["exp"]["dG"])
        }
    cc = CycleClosure.ccgraph_from_data(node_start, node_end, ddg, error, exp_data)
    cc.run()
    result_json["summary"]["pair"]["cc_ddG"] = []
    result_json["summary"]["pair"]["cc_error"] = []
    result_json["summary"]["node"]["dG"] = []
    result_json["summary"]["node"]["error"] = []

    for start, end in zip(node_start, node_end):
        key = f"{start}_to_{end}"
        edge = cc.graph.edges_dict[key]
        result_json["summary"]["pair"]["cc_ddG"].append(edge.get_data("cc_ddg"))
        result_json["summary"]["pair"]["cc_error"].append(edge.get_data("cc_ddg_error"))

    for key in result_json["summary"]["node"]["name"]:
        if key in node_start or key in node_end:
            node = cc.graph.nodes_dict[key]
            result_json["summary"]["node"]["dG"].append(node.get_data("cc_dg"))
            result_json["summary"]["node"]["error"].append(node.get_data("cc_dg_error"))
        else:
            result_json["summary"]["node"]["dG"].append(None)
            result_json["summary"]["node"]["error"].append(None)
    return json.dumps(result_json)


def run_cycle_closure(
    node_start: List[str],
    node_end: List[str],
    ddg: List[float],
    error: List[float],
    exp_data: Dict[str, float] = None,
    bias_node: List[str] = None,
) -> Tuple[Dict[str, Dict[str, float]], Dict[Any, Any]]:

    cc = CycleClosure.ccgraph_from_data(node_start, node_end, ddg, error, exp_data, bias_node)
    cc.run()
    result = defaultdict(dict)
    for edge, data in cc.graph.edges_iter(data=True):
        result[edge.name]["cc_ddg"] = data.get("cc_ddg")
        result[edge.name]["cc_ddg_error"] = data.get("cc_ddg_error")
    for node, data in cc.graph.nodes_iter(data=True):
        result[node.name]["cc_dg"] = data.get("cc_dg")
        result[node.name]["cc_dg_error"] = data.get("cc_dg_error")
    return result, cc.hystersis


@dataclass(order=True)
class PrioritizedItem:
    error: float
    dg: float
    path: List[str]

    def __iter__(self):
        return iter((self.error, self.dg, self.path))


def get_mid_node(graph):
    """Getting the middle node, which has the smallest distance to other nodes"""
    path_length = {k: {n: len(v) for n, v in n_v.items()} for k, n_v in nx.all_pairs_shortest_path(graph)}
    path_max_length = {k: max(v.values()) for k, v in path_length.items()}
    radius = min(path_max_length.values())
    centers = [n for n, r in path_max_length.items() if r == radius]
    centers_and_sum_dist = {x: sum(path_length[x].values()) for x in centers}
    min_sum_dist = min(centers_and_sum_dist.values())
    min_centers = [x for x in centers_and_sum_dist if centers_and_sum_dist[x] == min_sum_dist]
    return sorted(min_centers, key=lambda x: (-x.bias, x.name))[0]


def find_closed_paths(graph: Graph, cutoff=6):
    # For the cycle basis is not the smallest cycle, although we could use this to form other cycles
    # but we cannot use 3 <= edge_number <= 6 to save other cycles.
    # We could generated all cycles and filter next,
    # in this way, we cannot remove cycles during the molecule_dynamics
    # Maybe it is time consuming.
    # For this reason, we use DFS to find other cycles
    paths = []

    def find_node_circle(node):
        path = [node]
        stack = [graph.neighbors(node)]  # generator
        while stack:
            children = stack[-1]
            try:
                child = next(children)
                if len(path) < cutoff:
                    if child == node:
                        yield path + [node]
                    elif child not in path:
                        path.append(child)
                        stack.append(graph.neighbors(child))
                else:  # we now reach the cutoff, test children contains node or not
                    if child == node or node in children:  # children has consumed child
                        yield path + [node]
                    stack.pop()
                    path.pop()
            except StopIteration:  # finished
                stack.pop()
                path.pop()

    for node in graph.nodes_iter():
        for path in find_node_circle(node):
            if len(path) <= 3:
                continue
            paths.append(path[:-1])
    return paths


class CycleClosure:
    def __init__(self, graph=None, gpickle=None, ddg=None, dg=None):
        if graph:
            self.graph = graph
        else:
            if gpickle:
                self.graph: Graph = Graph.load(gpickle)
            if ddg:
                self.assign_ddg_data(ddg)
        self.hystersis = {}

    @classmethod
    def ccgraph_from_csv_or_dataframe(cls, csvfile_or_df, exp_dg_df=None, split_keyword="_to_"):
        cc = cls()
        graph = Graph()
        if isinstance(csvfile_or_df, str) or isinstance(csvfile_or_df, Path):
            df = pd.read_csv(csvfile_or_df)
        else:
            df = csvfile_or_df
        for i, edge in enumerate(df["name"]):
            left, right = edge.split(split_keyword)
            leftNode = Node(graph, left)
            rightNode = Node(graph, right)
            graph.add_node(leftNode)
            graph.add_node(rightNode)
            edge = graph.add_edge(leftNode, rightNode)
            edge.direction = (leftNode, rightNode)
            # make sure the rightNode and leftNode direction is correct
            graph[leftNode][rightNode]["ddg"] = df["ddg"][i]
            graph[leftNode][rightNode]["ddg_error"] = df["ddg_error"][i]

        if exp_dg_df is not None:
            for _, row in exp_dg_df.iterrows():
                if row["name"] in graph.nodes_dict:
                    graph.nodes_dict[row["name"]].set_data("dg", row["dg"])
        cc.graph = graph
        return cc

    @classmethod
    def ccgraph_from_data(cls, node_start, node_end, ddg, error, exp_data=None, bias_node=None):
        cc = cls()
        graph = Graph()
        for i in range(len(node_start)):
            leftNode = Node(graph, node_start[i])
            rightNode = Node(graph, node_end[i])
            graph.add_node(leftNode)
            graph.add_node(rightNode)
            edge = graph.add_edge(leftNode, rightNode)
            edge.direction = (leftNode, rightNode)
            graph[leftNode][rightNode]["ddg"] = ddg[i]
            graph[leftNode][rightNode]["ddg_error"] = error[i]
        cc.graph = graph
        nodes_set = set(node_start) | set(node_end)
        if exp_data is not None:
            for key, value in exp_data.items():
                if key in nodes_set:
                    cc.graph.nodes_dict[key].set_data("dg", value)
        if bias_node is not None:
            for node in bias_node:
                cc.graph.nodes_dict[node].bias = True
        return cc

    def find_node_cutoff(self):
        N = self.graph.number_of_nodes()
        if N > 20:
            possible_stars = []
            R = len([n for n in self.graph.nodes_iter() if n.bias])
            # N : number of nodes
            # R : number of reference nodes
            # D : node degree
            for node in self.graph.nodes_iter():
                D = node.degree
                if node.bias:  # bias node connect most other nodes
                    possible_stars.append(D >= max(R + (N - R) / 2, N - 1 - (N - R) / 2))
                else:  # other node connect few other nodes
                    possible_stars.append(((D >= R / 2) & D <= (min(R + (N - R) / 2, N - 1 - (N - R) / 2))))
            if all(possible_stars):
                self.cutoff = NODE_CUTOFF_STAR_GRAPH
                logger.info("As the graph is star-like, ")
                logger.info(f"Using the cutoff {self.cutoff} to calculate the error")
                return
        self.cutoff = NODE_CUTOFF
        logger.info(f"Using the cutoff {self.cutoff} to calculate the error")

    def assign_ddg_data(self, ddg):
        for edge in self.graph.edges_iter():
            if edge.name in ddg.keys():
                edge.set_data("ddg", ddg[edge.name]["ddg"])
                edge.set_data("ddg_error", ddg[edge.name]["ddg_error"])

    def generate_subgraph_from_ddg_data(self):
        self.ddg_graph = self.graph.edge_subgraph(
            [edge.nodes for edge in self.graph.edges_iter() if edge.get_data("ddg") is not None]
        )
        self.ddg_subgraphs = []
        for nodes in nx.connected_components(self.ddg_graph):
            subgraph = self.ddg_graph.subgraph(nodes)
            # if dg in graph, then using those dg to calcuate other dg
            # if dg not in graph, we find the mid node, and set the dg of this node to 0
            has_dg = False
            for node in subgraph.nodes():
                if subgraph.nodes[node].get("dg", None):
                    has_dg = True
                    break
            if not has_dg:
                mid_node = get_mid_node(subgraph)
                logger.info("As there is no reference (bias) node, choose {} as reference node".format(mid_node))
                subgraph.nodes[mid_node]["dg"] = 0
            self.ddg_subgraphs.append(subgraph)

    def run(self):
        logger.info("Running cycle closure for the whole graph ...")
        self.find_node_cutoff()
        self.generate_subgraph_from_ddg_data()
        for i, subgraph in enumerate(self.ddg_subgraphs):
            logger.info(f"processing {i+1} subgraph")
            subcc = SubgraphCycleClosure(subgraph, cutoff=self.cutoff)
            subcc.run()
            self.hystersis.update(subcc.hystersis)

    def pka_correction(self, pka_file):
        pka_correction = PkaTautomerCorrection(self.graph, pka_file=pka_file)
        corrected_dg = pka_correction.get_corrected_dg()
        for edge, data in self.graph.edges_iter(data=True):
            edge.set_data("raw_ddg", data["ddg"])
            edge.set_data("raw_ddg_error", data["ddg_error"])
            node0, node1 = edge.nodes_name
            data["ddg"] += corrected_dg[node1] - corrected_dg[node0]
        for node, data in self.graph.nodes_iter(data=True):
            node.set_data("raw_cc_dg", data["cc_dg"])
            node.set_data("raw_cc_dg_error", data["cc_dg_error"])
        self.run()


class SubgraphCycleClosure:
    def __init__(self, graph, exp_error=EXPERIMENT_ERROR, cutoff=NODE_CUTOFF):
        self.graph: Graph = graph
        self.num_edges = self.graph.number_of_edges()
        self.corrected_ddg = None
        self.cccm_error = [-1.0] * self.num_edges
        self.exp_error = exp_error
        self.cutoff = cutoff
        self.path_for_each_node = {}

    def calc_cccm(self, cycle_list):
        cycle_cccm = set()
        for cycle in cycle_list:
            row_idx = [0] * self.num_edges
            for i in range(len(cycle)):
                try:
                    idx = self.edge_index[(cycle[i], cycle[(i + 1) % len(cycle)])]
                    row_idx[idx] = 1
                except:
                    idx = self.edge_index[(cycle[(i + 1) % len(cycle)], cycle[i])]
                    row_idx[idx] = -1
            signed_cccm = tuple(row_idx)
            oppo_signed_cccm = tuple([-item for item in row_idx])
            if signed_cccm not in cycle_cccm and oppo_signed_cccm not in cycle_cccm:
                cycle_cccm.add(signed_cccm)
        return np.array([list(cccm) for cccm in cycle_cccm])

    def find_cycles_and_cal_cccm(self, cutoff=NODE_CUTOFF):
        # base cycle and cccm for graph
        logger.info("Finding cycles and calculating cccm ...")
        self.edge_index = {(u, v): i for i, (u, v) in enumerate(self.graph.edges_iter())}
        base_cycles_list = nx.cycle_basis(self.graph)
        all_cycles_list = find_closed_paths(self.graph, cutoff)
        self.base_cycles_cccm = self.calc_cccm(base_cycles_list)
        self.all_cycles_cccm = self.calc_cccm(all_cycles_list)

    def calc_ddg(self):
        # svd to find the right null space and then using linear least square to find ddg
        logger.info("Calculating ddg by svd ...")
        ybar = np.array([data["ddg"] for _, data in self.graph.edges_iter(data=True)])
        if self.base_cycles_cccm.size:
            _, d, VT = np.linalg.svd(self.base_cycles_cccm)
            rank = np.count_nonzero(d > TOLERANCE)
            null_space_basis = np.array(VT[rank:, :].T)
            yhat = np.dot(null_space_basis, np.dot(null_space_basis.T, ybar))
        else:
            yhat = ybar
        self.corrected_ddg = yhat
        self.cccm_error = np.abs(yhat - ybar)
        for i, edge in enumerate(self.graph.edges_iter()):
            edge.set_data("cc_ddg", yhat[i])

    def calc_ddg_error(self):
        logger.info("Calculating ddg error ...")
        self.hystersis = {}
        for row in self.all_cycles_cccm:
            edge_list = [edge for val, edge in zip(row, self.graph.edges_iter(data=False)) if val]
            ddg_list = [-val * data["ddg"] for val, (_, data) in zip(row, self.graph.edges_iter(data=True)) if val]
            e = abs(sum(ddg_list)) / np.sqrt(len(ddg_list))
            for idx, val in enumerate(row):
                if val != 0 and e > self.cccm_error[idx]:
                    self.cccm_error[idx] = e
            cycle_pairs = []
            for edge in edge_list:
                cycle_pairs.append((edge.nodes_name[0], edge.nodes_name[1]))
            # # for pair in cycle_pair
            begin_edge = cycle_pairs.pop()
            cycle_paths = [begin_edge]
            while len(cycle_pairs):
                for i in range(len(cycle_pairs)):
                    if begin_edge[0] in cycle_pairs[i] or begin_edge[1] in cycle_pairs[i]:
                        cycle_paths.append(cycle_pairs[i])
                        del cycle_pairs[i]
                        begin_edge = cycle_paths[-1]
                        break
            self.hystersis[frozenset(cycle_paths)] = e
        for idx, ((n0, n1), data) in enumerate(self.graph.edges_iter(data=True)):
            e = self.cccm_error[idx] = max(self.cccm_error[idx], data["ddg_error"])
            self.graph[n0][n1]["cc_ddg_error"] = e
            self.graph[n0][n1]["cc_ddg_error2"] = e * e

    def path_dg_error(self, paths, ddg_dict):
        dg, error = 0, 0
        for i in range(len(paths) - 1):
            dg += ddg_dict[(paths[i], paths[i + 1])]
            error += self.graph[paths[i]][paths[i + 1]]["cc_ddg_error2"]
        return dg, error

    def calc_dg_and_error(self):
        logger.info("Calculation dg and the dg error ...")
        n = self.graph.number_of_nodes()
        ddg_dict = {}
        for idx, (n0, n1) in enumerate(self.graph.edges_iter()):
            ddg_dict[(n0, n1)] = self.corrected_ddg[idx]
            ddg_dict[(n1, n0)] = -self.corrected_ddg[idx]
        reference_node = []
        for node, data in self.graph.nodes(data=True):
            if data.get("dg", None) is not None:
                reference_node.append(node)
        k = len(reference_node)
        pq = []
        if k != 1:
            paths = nx.shortest_path(self.graph, weight="cc_ddg_error2")
        for i in range(k):
            node_i = reference_node[i]
            dg, sq_error = 0, 0
            for j in range(k):
                node_j = reference_node[j]
                dg += self.graph.nodes[node_j]["dg"]
                sq_error += self.exp_error
                if i != j:
                    ddg, err = self.path_dg_error(paths[node_i][node_j], ddg_dict)
                    dg -= ddg
                    sq_error += err
            item = PrioritizedItem(sq_error / k**2, dg / k, [node_i])
            heapq.heappush(pq, item)

        # BFS for multiple source and multiple sink
        while n:
            error, dg, path = heapq.heappop(pq)
            node = path[-1]
            if node not in self.path_for_each_node:
                self.path_for_each_node[node] = {"error": np.sqrt(error), "dg": dg, "path": path}
                n -= 1
                for neighbor in nx.all_neighbors(self.graph, node):
                    if neighbor not in self.path_for_each_node:
                        item = PrioritizedItem(
                            error + self.graph[node][neighbor]["cc_ddg_error2"],
                            dg + ddg_dict[(node, neighbor)],
                            path + [neighbor],
                        )
                        heapq.heappush(pq, item)

    def run(self):
        self.find_cycles_and_cal_cccm(self.cutoff)
        self.calc_ddg()
        self.calc_ddg_error()
        self.calc_dg_and_error()
        self.generate_results()

    def generate_results(self):
        for node_key, value in self.path_for_each_node.items():
            node_key.set_data("cc_dg", value["dg"])
            path = [item.name for item in value["path"]]
            node_key.set_data("cc_dg_path", path)
            node_key.set_data("cc_dg_error", value["error"])

        # do not need cc_ddg_errro2 data
        for edge in self.graph.edges_iter():
            edge.del_data("cc_ddg_error2")


if __name__ == "__main__":
    file = sys.argv[1]
    exp_file = sys.argv[2]
    df = pd.read_csv(file)
    names = df["name"]
    start, end = [], []
    for name in names:
        token = name.split("-")
        start.append(token[0])
        end.append(token[1])
        # start.append("-".join(token[:2]))
        # end.append("-".join(token[2:]))
    ddg = df["rbfe"].tolist()
    exp_ddg = df["exp"].tolist()
    error = df["rbfe_std"].tolist()
    exp = {}
    with open(exp_file) as f:
        for line in f:
            if line.strip() == "":
                continue
            key, value = line.split()
            exp[key] = float(value)
    res, cc_hy = run_cycle_closure(start, end, ddg, error, exp)
    with open("dg.csv", "w") as f:
        f.write("name,calc_dg,exp_dg\n")
        for key, value in res.items():
            for subkey, subvalue in value.items():
                if subkey == "cc_dg":
                    f.write(f"{key},{subvalue},{exp[key]}\n")

    with open("ddg.csv", "w") as f:
        f.write("name,calc_ddg,exp_ddg\n")
        for name, sub_ddg, sub_exp_ddg in zip(names, ddg, exp_ddg):
            f.write(f"{name},{sub_ddg},{sub_exp_ddg}\n")
