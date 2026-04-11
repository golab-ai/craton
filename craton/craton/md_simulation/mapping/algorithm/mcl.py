import numpy as np
import scipy.sparse as sparse

from ..graph import Graph


class MCL:
    """
    Markov Cluster Algorithm, a fast and scalable unsupervised clustering
    algorithm for graphs.
    This algorithm is based on molecule_dynamics of stochastic flow in graphs.
    """

    def __init__(self, graph: Graph, max_iter: int = 200, inflation: int = 2, pruning_cutoff: float = 1e-4):

        self.nodes = list(graph.nodes())
        num_nodes = graph.number_of_nodes()
        name2index = {node.name: i for i, node in enumerate(self.nodes)}
        self.max_iter = max_iter
        self.inflation = inflation
        self.pruning_cutoff = pruning_cutoff / num_nodes

        associated_matrix = np.zeros((num_nodes, num_nodes), dtype="f8")
        for e in graph.edges_iter():
            n0, n1 = e.nodes
            i = name2index[n0.name]
            j = name2index[n1.name]
            associated_matrix[i, j] = associated_matrix[j, i] = e.similarity or 0

        np.fill_diagonal(associated_matrix, np.amax(associated_matrix))
        self.associated_matrix = associated_matrix

    # All `m`s are of `scipy.sparse.csc_matrix` type.
    # CSC stands for Compressed Sparse Column.
    def prune(self, m):
        data = m.data
        for i in range(len(data)):
            if data[i] < self.pruning_cutoff:
                data[i] = 0
        m.eliminate_zeros()
        return m

    def normalize(self, m):
        data = m.data
        indptr = m.indptr
        for i in range(len(indptr) - 1):
            s = 0
            for j in range(indptr[i], indptr[i + 1]):
                s += data[j]
            for j in range(indptr[i], indptr[i + 1]):
                data[j] /= s
        return m

    def inflate(self, m):
        data = m.data
        for i in range(len(data)):
            data[i] = pow(data[i], self.inflation)
        return self.normalize(m)

    def expand(self, m):
        return m.dot(m)

    def run(self):
        curr_matrix = self.normalize(sparse.csc_matrix(self.associated_matrix))
        for i in range(self.max_iter):
            prev_matrix = curr_matrix
            # `expand` returns a new matrix. So `prev_matrix` won't be mutated by
            # the "expand => inflate => prune" operations.
            curr_matrix = self.prune(self.inflate(self.expand(curr_matrix)))
            if (curr_matrix - prev_matrix).count_nonzero() == 0:
                break

        csr_matrix = curr_matrix.tocsr()
        indices = csr_matrix.indices
        indptr = csr_matrix.indptr
        clusters = set()
        for i in range(len(indptr) - 1):
            c = tuple(self.nodes[j] for j in indices[indptr[i] : indptr[i + 1]])
            if c:
                clusters.add(c)

        return list(clusters)
