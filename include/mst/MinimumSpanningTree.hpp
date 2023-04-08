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
    MinimumSpanningTree(NetworKit::Graph &graph);

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
 * The class for the Boruvka minimum spanning tree algorithm
 */
class BoruvkaMinimumSpanningTree final : public MinimumSpanningTree {
 public:
    using MinimumSpanningTree::MinimumSpanningTree;

    /**
     * Execute the Boruvka minimum spanning tree algorithm.
     */
    void run();
};

/**
 * @ingroup mst
 * The class for the Karger-Klein-Tarjan randomized minimum spanning tree algorithm
 */
class KargerKleinTarjanMinimumSpanningTree final : public MinimumSpanningTree {
 public:
    using MinimumSpanningTree::MinimumSpanningTree;

    /**
     * Execute the Karger-Klein-Tarjan randomized minimum spanning tree algorithm.
     */
    void run();
};

} /* namespace Koala */
