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
 * The base class for the minimum spanning tree algorithms.
 *
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
 * The class for the Kruskal minimum spanning tree algorithm
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
 * The class for the Prim minimum spanning tree algorithm
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
 * The class for the Boruvka minimum spanning tree algorithm
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
 * The class for the Karger-Klein-Tarjan randomized minimum spanning tree algorithm
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

class ChazelleRubinfeldTrevisanMinimumSpanningTree final : public MinimumSpanningTree {
 public:
    using MinimumSpanningTree::MinimumSpanningTree;
    /**
     * Execute Chazelle-Rubinfeld-Trevisan randomized minimum spanning tree weight algorithm.
     *
     * eps - a constant in (0, 0.5), used to bound running time and result accuracy.
     * The smaller the eps, the more accurate the output.
     *
     * w - maximum edge weight - the algorithm assumes all edges have weights from {1,...,w}
     * as the algorithm is sublinear, w cannot be determined at runtime.
     */
    void run(unsigned int w, float eps = 0.1);

    /**
     * Do not use this function, instead use run with parameters
     * This function is neccessary for the class to compile, but is not implemented.
     */
    void run();

    /**
     * The algorithm does not calculate minimum spanning tree, only it's approximate weight.
     * Thus getForest method throws an exception.
     */
    const NetworKit::Graph& getForest() const;

    /**
     * Get the approximate weight of minimum spanning tree, calculated in run() method.
    */
    NetworKit::edgeweight getTreeWeight() const;

 private:
    float calculate_approximate_degree(float eps) const;
    float calculate_approximate_ccs_count(
            float eps, NetworKit::count bfs_bound, unsigned int w_bound) const;
    NetworKit::edgeweight calculate_approximate_tree_weight(float eps, unsigned int w) const;
    NetworKit::edgeweight tree_weight = NetworKit::nullWeight;
};

class Chazelle2000MinimumSpanningTree final : public MinimumSpanningTree {
 public:
    using MinimumSpanningTree::MinimumSpanningTree;
    void run();

 private:
    NetworKit::Graph mst(NetworKit::Graph G, int t);
    NetworKit::Graph msf(NetworKit::Graph G, int t);
    std::tuple<NetworKit::Graph, EdgeMap, NetworKit::Graph>
        boruvka_steps(NetworKit::Graph G, int c);
    int vertices_on_level(int dz);

    // The paper bounds the number of bad edges by 8m'/c + d^3n' and requires
    // it to be at most m'/2 + d^3n'. Its smallest hierarchy target is S(t, 1)^3 = 8.
    static constexpr int C = 16;
    static constexpr int MIN_NUMBER_NODES = 8;
};

}  // namespace Koala
