/*
 * GreedyVertexColoring.hpp
 *
 *  Created on: 22.10.2021
 *      Author: Krzysztof Turowski
 */

#pragma once

#include <map>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace Koala {

/**
 * @ingroup coloring
 * The base class for the greedy vertex coloring heuristics.
 *
 */
class GreedyVertexColoring : public NetworKit::Algorithm {

public:
    /**
     * Given an input graph, set up the greedy vertex coloring procedure.
     *
     * @param graph The input graph.
     */
    GreedyVertexColoring(const NetworKit::Graph &graph);

    /**
     * Return the coloring found by the algorithm.
     *
     * @return a map from nodes to colors.
     */
    const std::map<NetworKit::node, int>& getColoring() const;

protected:
    const std::optional<NetworKit::Graph> graph;
    std::map<NetworKit::node, int> colors;

    std::map<NetworKit::node, int>::iterator greedy_color(NetworKit::node v);
};

/**
 * @ingroup coloring
 * The class for the random sequential greedy vertex coloring heuristic.
 */
class RandomSequentialVertexColoring final : public GreedyVertexColoring {

public:
    using GreedyVertexColoring::GreedyVertexColoring;

    /**
     * Execute the random sequential greedy vertex coloring heuristic.
     */
    void run();
};

/**
 * @ingroup coloring
 * The class for the largest first greedy vertex coloring heuristic.
 */
class LargestFirstVertexColoring final : public GreedyVertexColoring {

public:
    using GreedyVertexColoring::GreedyVertexColoring;

    /**
     * Execute the largest first greedy vertex coloring heuristic.
     */
    void run();

private:
    std::vector<NetworKit::node> largest_first_ordering();
};

} /* namespace Koala */
