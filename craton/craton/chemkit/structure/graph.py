import networkx as nx
from ...utils import logger

def make_mole_as_graph(molecule, ignore_elems=None):
    G = nx.Graph()
    if ignore_elems is None:
        ignore_elems = []
    atoms = [i for i in range(len(molecule.Atoms)) if molecule.Atoms[i].elem not in ignore_elems]
    G.add_nodes_from(atoms)
    G.add_edges_from([(bond.a1, bond.a2) for bond in molecule.Bonds if bond.a1 in atoms and bond.a2 in atoms])
    return G
def get_frag_graph(self, frag_type="sketch_frag", ignore_side=True):
    self.scaffoldG = nx.Graph()
    if not hasattr(self, frag_type):
        logger.warning("please run fragmentation before")
        return None
    frags = getattr(self, frag_type)
    edges = []
    self.scaffoldG.add_edges_from([index for index in frags])
    for index, frag in frags.items():
        if frag["type"] == "Side" and ignore_side:
            continue
        for rr in frag["frag_connects"]:
            if f"{index}-{rr[0]}" not in edges and f"{rr[0]}-{index}" not in edges:
                edges.append(f"{index}-{rr[0]}")
                label = "-".join(sorted([frag["type"], frags[rr[0]]["type"]]))
                label = frag["type"] + "-" + frags[rr[0]]["type"]
                self.scaffoldG.add_edge(index, rr[0], label=label)
def calc_bond_distance(self, p0, p1):
    return len(nx.shortest_path(self.G, source=p0, target=p1))


def in_same_smallest_ring(self, *ids):
    try:
        smallest_ring_sizes = [min(self.Atoms[i].ring_size) for i in ids]
    except ValueError:
        return False
    if len(set(smallest_ring_sizes)) > 1:
        return False
    for k, v in self.ring_dict.items():
        ring = v[:-1]
        if len(ring) == smallest_ring_sizes[0] and all(i in ring for i in ids):
            return True
    return False
