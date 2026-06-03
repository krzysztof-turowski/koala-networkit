/*
 * MinimumSpanningTree.hpp
 *
 *  Created on: 08.04.2023
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

#include <optional>
#include <tuple>
#include <unordered_map>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/structures/UnionFind.hpp>

namespace Koala {

using EdgeMap = std::unordered_map<NetworKit::Edge, NetworKit::Edge>;

/**
 * @ingroup mst
 * Base class for minimum spanning forest algorithms on undirected graphs.
 */
class MinimumSpanningTree : public NetworKit::Algorithm {
 public:
    /**
     * Given an input graph, set up the minimum spanning tree procedure.
     *
     * @param graph The input graph.
     */
    explicit MinimumSpanningTree(NetworKit::Graph &graph);

    /**
     * Return the spanning tree found by the algorithm.
     *
     * @return a spanning tree.
     */
    const NetworKit::Graph& getForest() const;

    /**
     * Verify the result found by the algorithm using O(n + m) MST verification algorithm
     * from Hagerup, An Even Simpler Linear-Time Algorithm for Verifying Minimum Spanning Trees.
     */
    void check() const;

 protected:
    NetworKit::Graph graph, tree;
    void initialize();
};

/**
 * @ingroup mst
 * Kruskal's greedy minimum spanning forest algorithm.
 *
 * The algorithm considers edges in nondecreasing order of weight and joins two components
 * whenever the selected edge does not create a cycle.
 *
 * @see https://doi.org/10.1090/S0002-9939-1956-0078686-7
 */
class KruskalMinimumSpanningTree final : public MinimumSpanningTree {
 public:
    using MinimumSpanningTree::MinimumSpanningTree;

    /**
     * Execute the Kruskal minimum spanning tree algorithm.
     */
    void run();
};

/**
 * @ingroup mst
 * Prim's greedy minimum spanning tree algorithm.
 *
 * Starting from one vertex, the algorithm repeatedly adds the cheapest edge leaving the
 * current tree.
 *
 * @see https://doi.org/10.1002/j.1538-7305.1957.tb01515.x
 */
class PrimMinimumSpanningTree final : public MinimumSpanningTree {
 public:
    using MinimumSpanningTree::MinimumSpanningTree;

    /**
     * Execute the Prim minimum spanning tree algorithm.
     */
    void run();
};

/**
 * @ingroup mst
 * Boruvka's component-contraction minimum spanning forest algorithm.
 *
 * Each phase chooses a cheapest outgoing edge for every current component and contracts the
 * resulting forest.
 */
class BoruvkaMinimumSpanningTree : public MinimumSpanningTree {
    friend class MinimumSpanningTree;
    friend class Chazelle2000MinimumSpanningTree;
 public:
    using MinimumSpanningTree::MinimumSpanningTree;

    /**
     * Execute the Boruvka minimum spanning tree algorithm.
     */
    void run();

 protected:
    static std::optional<NetworKit::Graph> iterate(
        NetworKit::Graph &G, NetworKit::Graph &F,
        NetworKit::UnionFind &union_find, EdgeMap &E,
        NetworKit::count steps, bool get_branching_tree);
};

/**
 * @ingroup mst
 * Randomized linear-time minimum spanning tree algorithm of Karger, Klein, and Tarjan.
 *
 * The algorithm combines Boruvka contractions, random sampling, and removal of edges that are
 * heavy with respect to a recursively computed forest.
 *
 * @see https://doi.org/10.1145/201019.201022
 */
class KargerKleinTarjanMinimumSpanningTree final : public BoruvkaMinimumSpanningTree {
 public:
    using BoruvkaMinimumSpanningTree::BoruvkaMinimumSpanningTree;

    /**
     * Execute the Karger-Klein-Tarjan randomized minimum spanning tree algorithm.
     */
    void run();

 protected:
    static void recurse(NetworKit::Graph &G, NetworKit::Graph &F);
    static void discard_random_edges(NetworKit::Graph &G, NetworKit::Graph &subgraph);
    static void remove_heavy_edges(NetworKit::Graph &G, NetworKit::Graph &subgraph);
};

/**
 * @ingroup mst
 * Sublinear-time MST-weight approximation algorithm of Chazelle, Rubinfeld, and Trevisan.
 *
 * For a connected adjacency-list graph with weights in {1, ..., w}, the algorithm estimates
 * the MST weight within relative error eps using randomized component-count estimates.
 *
 * @see https://doi.org/10.1137/S0097539702403244
 */
class ChazelleRubinfeldTrevisanMinimumSpanningTree final : public MinimumSpanningTree {
 public:
    using MinimumSpanningTree::MinimumSpanningTree;
    /**
     * Estimate the minimum spanning tree weight.
     *
     * @param w Maximum edge weight. The algorithm assumes weights from {1, ..., w}.
     * @param eps Relative-error bound in (0, 0.5).
     */
    void run(unsigned int w, float eps = 0.1);

    /**
     * Throw because the approximation algorithm requires explicit bounds.
     */
    void run();

    /**
     * Throw because the approximation algorithm computes only a weight estimate.
     */
    const NetworKit::Graph& getForest() const;

    /**
     * Return the approximate minimum spanning tree weight.
     *
     * @return The MST-weight estimate calculated by run().
     */
    NetworKit::edgeweight getTreeWeight() const;

 private:
    float calculate_approximate_degree(float eps) const;
    float calculate_approximate_ccs_count(
            float eps, NetworKit::count bfs_bound, unsigned int w_bound) const;
    NetworKit::edgeweight calculate_approximate_tree_weight(float eps, unsigned int w) const;
    NetworKit::edgeweight tree_weight = NetworKit::nullWeight;
};

/**
 * @ingroup mst
 * Chazelle's deterministic minimum spanning tree algorithm.
 *
 * The comparison-based algorithm uses Boruvka contractions and soft heaps to achieve
 * O(m alpha(m, n)) time, where alpha is an inverse-Ackermann function.
 *
 * @see https://doi.org/10.1145/355541.355562
 */
class Chazelle2000MinimumSpanningTree final : public MinimumSpanningTree {
 public:
    using MinimumSpanningTree::MinimumSpanningTree;

    /**
     * Execute Chazelle's deterministic minimum spanning tree algorithm.
     */
    void run();

 private:
    NetworKit::Graph mst(NetworKit::Graph G, int t);
    NetworKit::Graph msf(NetworKit::Graph G, int t);
    std::tuple<NetworKit::Graph, EdgeMap, NetworKit::Graph>
        boruvka_steps(NetworKit::Graph G, int c);
    int vertices_on_level(int dz);
};

}  // namespace Koala
