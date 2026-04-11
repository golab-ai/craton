import itertools
import random

import networkx


def _calc_pheromone(similog):
    return 1.5 - min(1.0, similog)


# def _to_cedge(edges):
#     # dummy to keep the same interface as ant.pyx
#     return edges


def _run_ant_colony(self, cutoff_rule, num_iter=100, max_best_hits=40):
    """
    Ant Colony algorithm to optimize the graph topology.

    :type  cutoff_rule: `similarity.Cutoff`
    :param cutoff_rule: The cutoff rule for similarity calculation
    :type     num_iter: `int`
    :param    num_iter: Number of iterations to analyze to optimize the graph.
    """
    nodes = self.nodes()
    edges = sorted(self.edges_iter(), key=lambda e: e.similarity)
    # This algorithm has an unfortunate dependency on the order of elements
    # in the list `edges', which makes the result hard to reproduce. The
    # culprit is using the `_gen_random_graph' function (see below). To
    # alleviate this problem, we here sort the `edges'. It doesn't matter
    # which key to use for sorting: As long as we have the same edge
    # sequence, we should get consistent optimization results.

    # Constants
    MAX_PATHLEN = 5
    NUM_ANTS = len(nodes)
    S1_TARGET = _connectivity(self, MAX_PATHLEN)

    # This algorithm doesn't use the similarity score directly.
    # Converts all similarity scores to log(similarity)/log(cutoff).
    for e in edges:
        setattr(e, "similog", cutoff_rule.similog(e.similarity))

    connected = {}
    for e in edges:
        if e.get_data("KEEP"):
            n0, n1 = e.nodes
            connected[n0] = connected.get(n0, 0) + 1
            connected[n1] = connected.get(n1, 0) + 1

    best_s = float("-inf")
    best_g = networkx.Graph()
    retries = 0
    target_num_edges = 1.5 * len(nodes)

    total_pheromone = sum([_calc_pheromone(e.similog) for e in edges])
    pheromones = []
    # Invariant: sum(pheromones) of the graph should always be 1.
    for e in edges:
        n0, n1 = e.nodes
        pheromones.append(
            _calc_pheromone(e.similog) / total_pheromone / (connected.get(n0, 0) + connected.get(n1, 0) + 1)
        )

    random.seed(0)
    while _connectivity(best_g, MAX_PATHLEN) < S1_TARGET and retries < 3:
        retries += 1
        best_hits = 0

        # Resets the pheromone values.
        for e, p in zip(edges, pheromones):
            setattr(e, "pheromone", p)

        for num in range(num_iter):
            # Generating random graphs is the slowest step.
            # For each generated graph, the `_gen_random_graph' function has
            # to evaluate the graph's connectivity, which is slow and scales
            # as N^3, where N is the number of nodes in the graph.
            score_ant = [_gen_random_graph(target_num_edges, edges, MAX_PATHLEN) for _ in range(NUM_ANTS)]

            scores = []
            for score, ant in score_ant:
                scores.append(score)
                if score > best_s:
                    best_g = ant
                    best_s = score
                    best_hits = 0
            best_hits += 1
            if best_hits > max_best_hits:
                # If we continuously hit the same best graph for 40 times,
                # we consider the optimization has converged.
                # Different numbers have been tested, 40 is the one that can
                # preserve the same result as before.
                break
            _tournament(score_ant)
            num_edges = 0
            for w, a in score_ant:
                n = a.number_of_edges()
                num_edges += w * n
                for _, __, data in a.edges(data=True):
                    e = data["orig_edge"]
                    e.pheromone += w / n
            for e in edges:
                e.pheromone /= 2

            # target_num_edges = old_div((target_num_edges + num_edges), 2)
            target_num_edges = (target_num_edges + num_edges) / 2

    self.del_edge_data("KEEP")
    for _, _, data in best_g.edges(data=True):
        e = data["orig_edge"]
        e.keep = True


def _count_simple_paths(g, n0, n1, max_pathlen):
    # FIXME: this has been replaced.
    # We should probably delete it sometime in the future.
    # Current usage is for comparing speed to cython code.

    path_count = 0
    count = 0
    # `all_simple_paths(...)' returns an iterator. This function scales as
    # O(N + E), where `N' and `E' are the number of nodes and edges,
    # respectively.
    for _ in networkx.all_simple_paths(g, n0, n1, max_pathlen):
        path_count += 1
        if path_count >= 2:
            break
    count += path_count
    return count


def _connectivity(g, max_pathlen=5):
    """
    For each node pair, we get the result of 1 if there is only 1 connecting
    path, or 2 if 2 or more connecting paths, or 0 if none. We return the sum of
    the results of all pairs.

    :type g: `networkx.Graph` or its subclass
    """
    # This whole function scales as O(N^3).
    count = 0
    for n0, n1 in itertools.combinations(g.nodes, 2):
        count += _count_simple_paths(g, n0, n1, max_pathlen)
    return count


def _gen_random_graph(target_num_edge, edges, max_pathlen):
    s = 0
    g = networkx.Graph()
    for edge in edges:
        keep = edge.get_data("KEEP")
        if keep or (random.random() < target_num_edge * edge.pheromone):
            g.add_edge(*edge.nodes, orig_edge=edge)
            s += edge.similog
    s1 = _connectivity(g, max_pathlen)
    return [2 * s1 - s, g]


def _tournament(score_ant):
    """
    Normalize the scores in the score-ant list.
    The normalization is done as follows:
    We convert a list of raw scores into a list of values that are evenly
    spaced, yet the relative rank of each raw score is kept the same by the
    normalized value. For example, say we have 5 ants with scores
    [4, 3, 5, 1, 2], the normalized values would be
    [7/25, 5/25, 9/25, 1/25, 3/25]. Note that sum of the normalized values
    in the list remains 1.

    :type  score_ant: `list` of lists. Each list element should have at least
                      one element, which should be of type `float`.
    :param score_ant: The first element of each list element is the raw score to
                      be normalized. This argument will be mutated: The raw
                      scores will be replaced by its normalized value.
    """
    n = len(score_ant)
    scale = 1.0 / (n * n)
    new_score_ant = list(enumerate(score_ant))  # elem = (index, [score, ant])
    new_score_ant.sort(key=lambda x: x[1][0])  # Sort based on score.
    for j, (i, _) in enumerate(new_score_ant):
        score_ant[i][0] = (2 * j + 1) * scale
