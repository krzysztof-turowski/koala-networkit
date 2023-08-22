/*
 * MinimumSpanningTree.hpp
 *
 *  Created on: 08.04.2023
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

#include <optional>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/structures/UnionFind.hpp>

namespace Koala {

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

 protected:
    std::optional<NetworKit::Graph> graph, tree;
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
 public:
    using MinimumSpanningTree::MinimumSpanningTree;

    /**
     * Execute the Boruvka minimum spanning tree algorithm.
     */
    void run();

 protected:
    static void iterate(
        NetworKit::Graph &G, NetworKit::Graph &F, NetworKit::UnionFind &union_find,
        NetworKit::count steps);
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

} /* namespace Koala */
