import itertools
import time

import matplotlib
import psutil
from rdkit import Chem

from ...utils import logger
#from compuchem.chemistry.format.rdkit_wrapper import RdKitWrapper
from .algorithm.ant_colony import _run_ant_colony
from .graph import Graph
from .algorithm.mcl import MCL
from .algorithm.similarity import Cutoff
#from compuchem.molecule_dynamics.algorithm.fep.util import create_ligand_from_sdf

matplotlib.use("Agg")

_MAX_ANT_NUM_ITERATIONS = 200
_DEFAULT_ANT_MAX_BEST_HITS = 40
_DEFAULT_NUM_PROCS = psutil.cpu_count(logical=False) - 1

class PairNetwork:
    def __init__(self,):
        pass

    @staticmethod
    def generate_bias_unbias_nodes(bias_nodes,gg):
        if isinstance(bias_nodes, str):
            bias_nodes = [bias_nodes]
        bias_nodes = list(set(bias_nodes) if bias_nodes else set([]))
        nodes_dict = gg.nodes_dict
        bias_nodes = [nodes_dict[node_name] for node_name in bias_nodes]

        for node in bias_nodes:
            node.bias = True
        return bias_nodes

    @staticmethod
    def make_full_graph(gg,nbunch):
        #logger.info("Create full pair graph ...")
        if nbunch:
            all_nodes = gg.nodes()
            all_node_pairs = set(itertools.combinations(all_nodes, 2))
            excluded_nodes = set(all_nodes) - set(nbunch)
            excluded_pairs = set(itertools.combinations(excluded_nodes, 2))
            node_pairs = all_node_pairs - excluded_pairs
        else:
            node_pairs = list(gg.node_pairs())
        gg.add_edges_from(node_pairs)

    @staticmethod
    def create_graph_from_pairs(gg,user_pair_list):
        logger.info("Create graph from pairs ...")
        name_node_dict = gg.nodes_dict
        pairs = [(name_node_dict[pair[0]], name_node_dict[pair[1]]) for pair in user_pair_list]
        gg.add_edges_from(pairs, direction=True)

    @staticmethod
    def make_star_graph(gg,bias_nodes):
        logger.info(f"Create {bias_nodes} star graph ...")
        # nbunch = self.nbunch or self.graph.nodes()
        for node in gg.nodes_iter():
            for node_b in bias_nodes:
                #if not gg.has_edge(node, node_b):
                if not gg.has_edge([node,node_b]):
                    edge = gg.add_edge(node, node_b)
                    if edge:
                        edge.direction = (node_b, node)

    @staticmethod
    def create_graph_from_molecules(ligands,
                                    topology="normal",
                                    user_pair_list=None,
                                    bias_nodes=None,
                                    core=None,
                                    nbunch=None):
        gg = Graph()
        gg.add_nodes_from(ligands)
        bias_nodes = PairNetwork.generate_bias_unbias_nodes(bias_nodes,gg)
        gg.bias_nodes = bias_nodes
        if topology in ["normal","full"]:
            PairNetwork.make_full_graph(gg,nbunch)
        elif topology == "star":
            PairNetwork.make_star_graph(gg,bias_nodes)
        elif topology == "user_pairs":
            PairNetwork.create_graph_from_pairs(gg,user_pair_list)
        return gg

    @staticmethod
    def generate_core_atoms(gg,core):
        if core:
            logger.info(f"Select core atoms using the smile pattern {core}")
            core = Chem.MolFromSmiles(core)
            core = Chem.AddHs(core, explicitOnly=True)
            for name, node in gg.nodes_dict.items():
                rd_mol = RdKitWrapper(node.struct).molecule
                if core_atoms := rd_mol.GetSubstructMatch(core):
                    node.core = set(core_atoms)
                else:
                    logger.warning(f"{name} cannot match the input smiles, using all atoms instead !")
                    node.core = set(range(len(node.struct.Atoms)))        
    
    @staticmethod
    def reduce_normal_graph(gg,nbunch=None,bias_nodes=None,):
        nbunch = nbunch or gg.nodes
        for node in gg.nodes_iter():
            max_heap = []
            for node_b in gg.bias_nodes:
                if node_b is not node:
                    max_heap.append((gg.get_edge(node, node_b).similarity, node_b))
            if max_heap:
                max_heap.sort(reverse=True)
                start_node = node
                hub_node = max_heap[0][1]
                edge = gg.get_edge(start_node, hub_node)
                edge.keep = True
                edge.direction = (start_node, hub_node)

        gg.remove_edges_from(e for e in gg.edges_iter() if e.similarity == 0)
        clusters = MCL(gg).run()

        #logger.info(f"Clusters in Graph: {len(clusters)}")
        cutoffRule = Cutoff(cutoff=0)

        start_time = time.time()
        #logger.info("Start ant colony optimization algorithm...")
        for i, e in enumerate(clusters, start=1):
            #logger.info("Optimizing the subgraph of cluster#%d..." % i)
            #logger.info(f"{len(e)} nodes")
            subgraph = gg.subgraph(e)
            subgraph_size = len(subgraph)
            if subgraph_size > 2:
                _run_ant_colony(
                    subgraph, cutoffRule, num_iter=_MAX_ANT_NUM_ITERATIONS, max_best_hits=_DEFAULT_ANT_MAX_BEST_HITS
                )
            elif subgraph_size == 2:
                the_only_edge = next(subgraph.edges_iter())
                the_only_edge.set_data("KEEP", True)
            #logger.info("Optimizing the subgraph of cluster#%d... Done" % i)

        end_time = time.time()
        #logger.info(
        #    "Ant colony algorithm done (%.2f seconds). %d edges" % (end_time - start_time, gg.number_of_edges())
        #)

        gg.connect_clusters(clusters)
        gg.remove_bad_edges(check_tag=True)
        gg.del_edge_data("KEEP")


    @staticmethod
    def calculate_atom_mapping(gg):
        gg.calc_atom_mapping()

    @staticmethod
    def calculate_similarity(gg):
        gg.calc_similarity()

    @staticmethod
    def report_graph_attributes(gg,topology="normal"):
        simi_scores = [edge.similarity for edge in gg.edges_iter()]
        average_score = sum(simi_scores) / len(simi_scores)
        logger.debug(f"simi_scores: {simi_scores}")
        logger.info(f"total edges of the final graph: {len(simi_scores)}")
        logger.info(f"average similarity score of the final graph: {average_score:.3f}")
        gg.graph["average_similarity"] = average_score
        gg.graph["graph_topology"] = topology


class GraphGenerator_old:
    def __init__(self, nbunch=None, topology="normal", bias_nodes=None, num_procs=None, core=None):
        self.topology = topology
        self.bias_nodes = bias_nodes
        self.num_procs = num_procs or _DEFAULT_NUM_PROCS
        self.nbunch = nbunch
        self.core = core
        self.graph: Graph = None

    @classmethod
    def create_graph_from_gpickle(cls, gpickle_file):
        logger.info("Reading graph from gpickle file ...")
        gg = cls()
        gg.graph: Graph = Graph.load(gpickle_file)
        gg.report_graph_attributes()
        return gg

    @classmethod
    def create_graph_from_sdffile(
        cls, sdffile, topology="normal", nbunch=None, num_procs=None, bias_nodes=None, core=None
    ):
        num_procs = num_procs or _DEFAULT_NUM_PROCS
        logger.info("Start creating graph ...")
        gg = cls(nbunch, topology, bias_nodes, num_procs, core)
        inf_mols = create_ligand_from_sdf(sdffile)
        gg.graph: Graph = Graph()
        gg.graph.add_nodes_from(inf_mols)
        gg.generate_bias_unbias_nodes()
        gg.generate_core_atoms()

        if topology == "normal":
            gg.make_normal_graph()
        elif topology == "star":
            gg.make_star_graph()
        elif topology == "full":
            gg.make_full_graph()
        else:
            raise ValueError("The topol of graph must be 'normal', 'star', and 'full'")
        gg.report_graph_attributes()
        return gg

    @classmethod
    def create_graph_from_mols(cls, mols, topology="normal", nbunch=None, num_procs=1, bias_nodes=None, core=None):
        logger.info("Start creating graph ...")
        gg = cls(nbunch, topology, bias_nodes, num_procs, core)
        gg.graph: Graph = Graph()
        gg.graph.add_nodes_from(mols)
        gg.generate_bias_unbias_nodes()
        gg.generate_core_atoms()

        if topology == "normal":
            gg.make_normal_graph()
        elif topology == "star":
            gg.make_star_graph()
        elif topology == "full":
            gg.make_full_graph()
        else:
            raise ValueError("The topol of graph must be 'normal', 'star', and 'full'")
        gg.report_graph_attributes()
        return gg

    @classmethod
    def create_graph_from_pairs(
        cls, pairs, mols, topology="defined", num_procs=None, bias_nodes=None, core=None
    ) -> Graph:
        num_procs = num_procs or _DEFAULT_NUM_PROCS
        logger.info("Create graph from pairs ...")
        gg = cls(topology=topology, num_procs=num_procs, bias_nodes=bias_nodes, core=core)

        gg.graph: Graph = Graph()
        gg.graph.add_nodes_from(mols)
        name_node_dict = gg.graph.nodes_dict
        gg.generate_bias_unbias_nodes()
        gg.generate_core_atoms()

        pairs = [(name_node_dict[pair[0]], name_node_dict[pair[1]]) for pair in pairs]
        gg.graph.add_edges_from(pairs, direction=True)
        gg.graph.calc_atom_mapping(num_procs=gg.num_procs)
        gg.graph.calc_similarity(num_procs=gg.num_procs)
        gg.report_graph_attributes()
        return gg

    def make_full_graph(self):
        if self.nbunch:
            all_nodes = self.graph.nodes()
            all_node_pairs = set(itertools.combinations(all_nodes, 2))
            excluded_nodes = set(all_nodes) - set(self.nbunch)
            excluded_pairs = set(itertools.combinations(excluded_nodes, 2))
            node_pairs = all_node_pairs - excluded_pairs
        else:
            node_pairs = list(self.graph.node_pairs())
        self.graph.add_edges_from(node_pairs)
        self.graph.calc_atom_mapping(num_procs=self.num_procs)
        self.graph.calc_similarity(num_procs=self.num_procs)

    def generate_bias_unbias_nodes(self):
        if isinstance(self.bias_nodes, str):
            self.bias_nodes = [self.bias_nodes]
        self.bias_nodes = list(set(self.bias_nodes) if self.bias_nodes else set([]))
        nodes_dict = self.graph.nodes_dict
        self.bias_nodes = [nodes_dict[node_name] for node_name in self.bias_nodes]

        for node in self.bias_nodes:
            node.bias = True

    def generate_core_atoms(self):
        if self.core:
            logger.info(f"Select core atoms using the smile pattern {self.core}")
            core = Chem.MolFromSmiles(self.core)
            core = Chem.AddHs(core, explicitOnly=True)
            for name, node in self.graph.nodes_dict.items():
                rd_mol = RdKitWrapper(node.struct).molecule
                if core_atoms := rd_mol.GetSubstructMatch(core):
                    node.core = set(core_atoms)
                else:
                    logger.warning(f"{name} cannot match the input smiles, using all atoms instead !")
                    node.core = set(range(len(node.struct.Atoms)))

    def make_star_graph(self):
        # nbunch = self.nbunch or self.graph.nodes()
        for node in self.graph.nodes_iter():
            for node_b in self.bias_nodes:
                if not self.graph.has_edge(node, node_b):
                    edge = self.graph.add_edge(node, node_b)
                    if edge:
                        edge.direction = (node_b, node)
        self.graph.calc_atom_mapping(num_procs=self.num_procs)
        self.graph.calc_similarity(num_procs=self.num_procs)

    def make_normal_graph(self):
        self.make_full_graph()
        nbunch = self.nbunch or self.graph.nodes

        for node in self.graph.nodes_iter():
            max_heap = []
            for node_b in self.bias_nodes:
                if node_b is not node:
                    max_heap.append((self.graph.get_edge(node, node_b).similarity, node_b))
            if max_heap:
                max_heap.sort(reverse=True)
                start_node = node
                hub_node = max_heap[0][1]
                edge = self.graph.get_edge(start_node, hub_node)
                edge.keep = True
                edge.direction = (start_node, hub_node)

        self.graph.remove_edges_from(e for e in self.graph.edges_iter() if e.similarity == 0)
        clusters = MCL(self.graph).run()

        logger.info(f"Clusters in Graph: {len(clusters)}")
        cutoffRule = Cutoff(cutoff=0)

        start_time = time.time()
        logger.info("Start ant colony optimization algorithm...")
        for i, e in enumerate(clusters, start=1):
            logger.info("Optimizing the subgraph of cluster#%d..." % i)
            logger.info(f"{len(e)} nodes")
            subgraph = self.graph.subgraph(e)
            subgraph_size = len(subgraph)
            if subgraph_size > 2:
                _run_ant_colony(
                    subgraph, cutoffRule, num_iter=_MAX_ANT_NUM_ITERATIONS, max_best_hits=_DEFAULT_ANT_MAX_BEST_HITS
                )
            elif subgraph_size == 2:
                the_only_edge = next(subgraph.edges_iter())
                the_only_edge.set_data("KEEP", True)
            logger.info("Optimizing the subgraph of cluster#%d... Done" % i)

        end_time = time.time()
        logger.info(
            "Ant colony algorithm done (%.2f seconds). %d edges" % (end_time - start_time, self.graph.number_of_edges())
        )

        self.graph.connect_clusters(clusters)
        self.graph.remove_bad_edges(check_tag=True)
        self.graph.del_edge_data("KEEP")

    def report_graph_attributes(self):
        simi_scores = [edge.similarity for edge in self.graph.edges_iter()]
        average_score = sum(simi_scores) / len(simi_scores)
        logger.debug(f"simi_scores: {simi_scores}")
        logger.info(f"total edges of the final graph: {len(simi_scores)}")
        logger.info(f"average similarity score of the final graph: {average_score:.3f}")
        self.graph.graph["average_similarity"] = average_score
        self.graph.graph["graph_topology"] = self.topology

    @staticmethod
    def pairnetwork_show(gg,filename=None,edge_label="similarity"):
        gg.show_network(filename=filename, edge_label=edge_label)

    @staticmethod
    def atom_mapping_show(gg):
        for edge in list(gg.edges_iter()):
            matcher = edge.atom_mapping