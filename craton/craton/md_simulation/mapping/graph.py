import copy
import itertools
import weakref
from typing import Any, Dict, Generator, Tuple, Union

import matplotlib.pyplot as plt
import networkx

from ...utils import logger
#from .algorithm.dual_topology import assign_dual_topology
from .algorithm.mapping import assign_atom_mapping
from .algorithm.similarity import assign_similarity


def _property(key, doc=None):
    return property(
        lambda self: self._get(key), lambda self, value: self._set(key, value), lambda self: self._del(key), doc=doc
    )


def _overload(func):
    def wrapper(self, *arg, **kwargs):
        parent_obj = super(self.__class__, self)
        parent_func = getattr(parent_obj, func.__name__)
        try:
            return func(self, *arg, **kwargs)
        except TypeError as e:
            if str(e).startswith("%s()" % func.__name__):
                try:
                    return parent_func(*arg, **kwargs)
                except TypeError as e:
                    if str(e).startswith("%s() takes" % func.__name__):
                        raise RuntimeError(
                            "overloading resolution failed:\n"
                            "function name: %s\n"
                            "arg: %s\n"
                            "kwarg: %s\n" % (func.__name__, arg, kwargs)
                        )
            raise e

    return wrapper


class EdgeNameDict(dict):
    def __missing__(self, key):
        node1, node2 = key.split("_to_")
        edge_name = "_to_".join([node2, node1])
        return self[edge_name]


class Node(object):
    def __init__(self, graph, name):
        self._name = name
        self._graph = weakref.ref(graph)

    def __eq__(self, rhs):
        try:
            return self._name == rhs._name and self.graph is rhs.graph
        except Exception:
            return False

    def __lt__(self, rhs):
        return self._name < rhs._name

    def __repr__(self):
        return f"Node: {self._name}"

    def __hash__(self):
        return hash(self._name)

    def _set(self, key, value):
        self.graph._set_node_data(self, key, value)

    def _del(self, key):
        self.graph.del_node_data(self, key)

    def _get(self, key, default_value=None):
        graph = self.graph
        if graph is None:
            return default_value
        return graph._get_node_data(self, key, default_value)

    def get_data(self, key, default_value=None):
        return self._get(key, default_value)

    def set_data(self, key, value):
        self._set(key, value)

    @property
    def degree(self):
        return self.graph.degree(self)

    @property
    def name(self):
        return self._name

    @property
    def graph(self):
        return self._graph()

    @property
    def neighbors(self):
        return self.graph.neighbors(self)

    @property
    def struct(self):
        return self._get("STRUCTURE")

    @property
    def bias(self):
        return self._get("BIAS", False)

    @bias.setter
    def bias(self, value=None):
        if value is None:
            self._set("BIAS", False)
        else:
            self._set("BIAS", True)

    @property
    def core(self):
        return self._get("CORE")

    @core.setter
    def core(self, value=None):
        self._set("CORE", value)


class Edge(object):
    def __init__(self, graph, n0, n1):
        self._nodes = (n1, n0) if n0.name > n1.name else (n0, n1)
        self._graph = weakref.ref(graph)
        self._cached_hash = None
        self._is_core_hopping = False

    def __iter__(self):
        return iter(self.nodes)

    def __getitem__(self, i):
        return self._nodes[i]

    def _get(self, key, default_value=None):
        return self.graph._get_edge_data(self._nodes[0], self._nodes[1], key, default_value)

    def __hash__(self):
        if not self._cached_hash:
            self._cached_hash = hash(self.name)
        return self._cached_hash

    def _set(self, key, value):
        self.graph._set_edge_data(self._nodes[0], self._nodes[1], key, value)

    def _del(self, key):
        self.graph._del_edge_data(self._nodes[0], self._nodes[1], key)

    def get_data(self, key, default_value=None):
        return self._get(key, default_value)

    def set_data(self, key, value):
        self._set(key, value)

    def del_data(self, key):
        self._del(key)

    @property
    def graph(self):
        return self._graph()

    @property
    def nodes(self):
        direction = self._get("DIRECTION")
        if direction and direction[0] == self._nodes[1].name:
            return self._nodes[1], self._nodes[0]
        return self._nodes

    @property
    def structs(self):
        n0, n1 = self.nodes
        return n0.struct, n1.struct

    @property
    def name(self):
        node0, node1 = self.nodes
        return f"{node0.name}_to_{node1.name}"

    @property
    def nodes_name(self):
        node0, node1 = self.nodes
        return (node0.name, node1.name)

    @property
    def direction(self):
        return self.nodes[0].name, self.nodes[1].name

    @direction.setter
    def direction(self, node_pair):
        from_node, to_node = node_pair
        self._set("DIRECTION", (from_node.name, to_node.name))

    atom_mapping = _property("ATOM_MAPPING")
    similarity = _property("SIMILARITY")
    similarity_detail = _property("SIMILARITY_DETAIL")
    is_core_hopping = _property("IS_CORE_HOPPING")
    is_charge_hopping = _property("IS_CHARGE_HOPPING")
    dual_atom_mapping = _property("DUAL_ATOM_MAPPING")  # Dual Topology Atom Mapping
    keep = _property("KEEP")
    bad_edge = _property("BAD_EDGE")
    rest2_ligand_hot_atoms = _property("HOT_ATOM")
    rest2_ligand_hot_torsions = _property("HOT_TORSION")
    rest2_protein_hot_atoms = _property("PROTEIN_HOT_ATOMS")
    rest2_protein_hot_torsions = _property("PROTEIN_HOT_TORSIONS")
    rest2_scale_factor = _property("SCALE_FACTOR")

    @classmethod
    def validate(cls, n0, n1):
        """
        Check if the two given nodes can form a valid edge. Validity criteria:
          1. Atom mapping of the two structures must not be empty.
        (more might be added in the future).

        Note: Would be slow if overusing it.

        Example:

            n0, n1 = g.add_nodes_from([ct0, ct1])
            if Edge.validate(n0, n1):
                g.add_edge(n0, n1)
            else:
                print("Empty atom mapping")
        """
        result, msg = True, ""
        g = n0.graph
        edge = Edge(g, n0, n1)
        added_new_edge = False

        if not g.has_edge(edge):
            added_new_edge = True
            edge = g.add_edge(n0, n1)
            # edge.annotate()

        atom_mapping = edge.atom_mapping or {}

        if edge.similarity == 0.0:
            reason = ", ".join(rule for rule, sim in edge.similarity_chain if sim < 1e-2)
            result, msg = False, f"Edge are not similiar (Caused by {reason} rule(s))."

        if not all(atom_mapping.values()):
            result, msg = False, "no mapping data"

        if added_new_edge:
            g.remove_edge(edge)
        return result, msg


class Graph(networkx.Graph):
    def __init__(self, *args, **kwargs):
        super().__init__(self, *args, **kwargs)
        self.nodes_name_dict = dict()
        self.edges_name_dict = EdgeNameDict()

    def copy(self):
        return self.__copy__()

    def __copy__(self):
        new_graph = Graph()
        for node, data in self.nodes(data=True):
            new_node = Node(new_graph, node.name)
            new_graph.add_node(new_node)
            new_graph._node[new_node] = data
        for edge, data in self.edges_iter(data=True):
            old_n0, old_n1 = edge.nodes
            n0 = Node(new_graph, old_n0.name)
            n1 = Node(new_graph, old_n1.name)
            new_graph.add_edge(n0, n1)
            new_graph._adj[n0][n1] = data
            new_graph._adj[n1][n0] = data
        return new_graph

    def _get_node_data(self, node, key, default_value):
        return self.nodes[node].get(key, default_value)

    def _set_node_data(self, node, key, value):
        self.nodes[node][key] = value

    def _del_node_data(self, node, key):
        del self.nodes[node][key]

    def _get_edge_data(self, n0, n1, key, default_value=None):
        try:
            e = self[n0][n1]
        except KeyError:
            if n0.graph is n1.graph:
                raise KeyError("No edge between %s and %s" % (n0.name, n1.name))
            raise KeyError("No edge. Two nodes (%s, %s) are actaully from " "different graphs." % (n0.name, n1.name))

        return e.get(key, default_value)

    def _set_edge_data(self, n0, n1, key, value):
        self[n0][n1][key] = value

    def _del_edge_data(self, n0, n1, key):
        del self[n0][n1][key]

    # add_node, remove_node, add_edge, remove_edge
    def add_nodes_from(self, mols):
        return list(map(self.add_node_from, mols))

    def add_node_from(self, mol):
        new_node = Node(self, mol.mole_name)
        attr = {"STRUCTURE": mol}
        if hasattr(mol, "dg"):
            attr["dg"] = mol.dg
        networkx.Graph.add_node(self, new_node, **attr)
        return new_node

    def add_edge(self, n0, n1, **attr):
        if n0 not in self or n1 not in self:
            raise ValueError("Node(s) does not exist")
        if n0.name == n1.name:
            logger.warning(f"Trying building edge with same node {n0.name}, skip !")
            return
        networkx.Graph.add_edge(self, n0, n1, **attr)
        e = Edge(self, n0, n1)
        e.direction = e.nodes
        return e

    def add_edges_from(self, ebunch, direction=False, **attr):
        for e in ebunch:
            assert e[0] in self
            assert e[1] in self
        networkx.Graph.add_edges_from(self, ebunch, **attr)
        edges = [Edge(self, *e) for e in ebunch]
        for e_add, e in zip(edges, ebunch):
            if direction:
                if e_add.nodes[0] == e[0]:
                    e_add.direction = e_add.nodes
                else:
                    e_add.direction = e[0], e[1]
            else:
                e_add.direction = e_add.nodes
        return edges

    def get_edge(self, n0, n1):
        return Edge(self, n0, n1)

    def add_edge_data(self, key, value, nbunch=None):
        for _, _, data in super().edges(nbunch, data=True):
            data[key] = value

    def del_edge_data(self, key, nbunch=None):
        for _, _, data in super().edges(nbunch, data=True):
            data.pop(key, None)

    def node_pairs(self):
        return itertools.combinations(self.nodes(), 2)

    def edges_iter(
        self, nbunch=None, data=False, default=None
    ) -> Union[Generator[Edge, None, None], Generator[Tuple[Edge, Any], None, None]]:
        edges_iter = super(Graph, self).edges
        if data:
            ###for e in edges_iter(nbunch, data, default):
            for e in edges_iter(nbunch, data):
                yield Edge(self, e[0], e[1]), e[2]
        else:
            ###for e in edges_iter(nbunch, data, default):
            for e in edges_iter(nbunch, data):
                yield Edge(self, *e)

    def nodes_iter(self, *args, **kwargs):
        return iter(self.nodes(*args, **kwargs))

    @property
    def nodes_dict(self) -> Dict[str, Node]:
        if len(self.nodes_name_dict) == 0:
            self.nodes_name_dict = {node.name: node for node in self.nodes_iter()}
        return self.nodes_name_dict

    @property
    def edges_dict(self):
        if len(self.edges_name_dict) == 0:
            self.edges_name_dict = EdgeNameDict({edge.name: edge for edge in self.edges_iter()})
        return self.edges_name_dict

    @_overload
    def has_edge(self, e):
        try:
            return networkx.Graph.has_edge(self, *e)
        except KeyError:
            return False

    @_overload
    def remove_edge(self, e):
        super().remove_edge(*e)

    def calc_atom_mapping(self, num_procs=1, core=None):
        # The way of schrodigner do atom mapping is a bit complicated
        # Our implemention is much simple
        # 1. we have no core smarts
        # 2. we do not find cores, so the cores is the all of atoms in the structure
        # 3. As we have no core smarts, so the chirality check is not preformed (by count the atom number and core idx in core is equal)
        # 4. When we have core smarts, we first analyze the chirality check to return a atom mapping data, and do the
        # atom mapping for the next atoms.
        # TODO: SPECIFY core smarts in our program, at now, we skip this.
        # Core smiles is suppported now, but only one core
        
        assign_atom_mapping(list(self.edges_iter()))

    def calc_similarity(self, num_procs=1):
        # We omit calculate some rule
        # 1. softBond (Finished)
        # 2. SnapCoreRMSD, instead we use RMSD
        # 3. BidirectionSnapCore
        assign_similarity(list(self.edges_iter()), num_procs=num_procs)

    #def calc_dual_topology(self, num_procs=8):
        # get dual topology for the edge
    #    assign_dual_topology(list(self.edges_iter()), num_procs=num_procs)

    def write(self, filename):
        import pickle
        # write key result into a datafile
        g = networkx.Graph()
        g.graph = self.graph
        for node, data in self.nodes(data=True):
            g.add_node(node.name, **data)
        for n0, n1, data in self.edges(data=True):
            g.add_edge(n0.name, n1.name, **data)
        with open(filename, 'wb') as f:
            pickle.dump(g, f, pickle.HIGHEST_PROTOCOL)

    @classmethod
    def load(cls, pickle_file):
        import pickle
        with open(pickle_file, 'rb') as f:
            g = pickle.load(f)
        graph = cls()
        graph.graph = g.graph
        for node, data in g.nodes(data=True):
            new_node = Node(graph, node)
            graph.add_node(new_node, **data)
        for n0, n1, data in g.edges(data=True):
            for node in graph.nodes:
                if node.name == n0:
                    n0_node = node
                if node.name == n1:
                    n1_node = node
            graph.add_edge(n0_node, n1_node, **data)
        return graph

    def show_network(self, filename=None, edge_label="similarity"):
        labels = {}
        for node in self.nodes:
            labels[node] = node.name
        pos = networkx.spring_layout(self)
        networkx.draw(self, with_labels=True, labels=labels, pos=pos)
        # edge_labels = networkx.get_edge_attributes(self, 'similarity')
        edge_labels = {}
        for edge in self.edges_iter():
            if edge_label == "similarity":
                edge_labels[edge.nodes] = f"{edge.name} \n {edge.similarity:.3f}"
                networkx.draw_networkx_edge_labels(self, pos, edge_labels=edge_labels)
            elif edge_label == "ddg":
                edge_labels[edge.nodes] = f"{edge.name} \n {edge.get_data('ddg', 0):.3f}"
                edge_labels[edge.nodes] = f"{edge.name} \n {edge.get_data('ddg', 0):.3f}"

        if filename:
            plt.savefig(filename)
        else:
            plt.savefig("network.png")
        plt.show()

    def connect_clusters(self, clusters):
        """
        We use edges with the highest similarity scores to connect clusters.

        :type   clusters: `list` of `list`s of `Node` objects
        :param  clusters: List of clusters, where each one is a list of nodes
                          that are connected to each other within the subgraph.
                          No edges between clusters.

        :return: Number of clusters after the connection.
        """
        num_clusters = len(clusters)
        if num_clusters == 1:
            return num_clusters

        for i in range(num_clusters):
            this = clusters[i]
            others = copy.copy(clusters)
            others.pop(i)
            if len(this) == 2:
                for node in this:
                    self._connect_clusters_impl([node], i, others, 1)
            else:
                self._connect_clusters_impl(this, i, others)

        # Reevaluates the clusters.
        deleted_edge_data = {}
        for e, data in list(self.edges_iter(data=True)):
            # if "KEEP" not in data:
            if not e.keep:
                deleted_edge_data[e.nodes] = data
                self.remove_edge(e)
        clusters = list(networkx.connected_components(self))
        for (n0, n1), data in deleted_edge_data.items():
            self.add_edge(n0, n1, **data)

        if num_clusters == len(clusters):
            return num_clusters
        return self.connect_clusters(clusters)

    def _connect_clusters_impl(self, this, i, others, count=2):
        def simi(node_pair):
            return self.get_edge(*node_pair).similarity

        c2c_edges = []
        for e in others:
            c2c_edges.extend(networkx.edge_boundary(self, this, nbunch2=e))

        if len(c2c_edges) == 0:
            logger.warning("WARNING: Cannot connect cluster#%d with the others." % i)
            return

        c2c_edges = sorted(c2c_edges, key=simi, reverse=True)
        for node_pair in c2c_edges[:count]:
            self.get_edge(*node_pair).set_data("KEEP", True)

    def remove_bad_edges(self, ebunch=None, check_tag=False):
        """Remove any existing edge if it does not have valid atom mapping
        or it has 0.0 similarity.
        """
        edges = ebunch or list(self.edges_iter())
        for e in edges:
            if (check_tag and not e.get_data("KEEP")) or not e.validate(*e)[0]:
                self.remove_edge(e)
