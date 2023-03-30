/*
 * VertexColoring.hpp
 *
 *  Created on: 30.03.2023
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

#include <map>
#include <optional>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace Koala {

/**
 * @ingroup coloring
 * The base class for the vertex coloring algorithms.
 *
 */
class VertexColoring : public NetworKit::Algorithm {

public:
    /**
     * Given an input graph, set up the vertex coloring procedure.
     *
     * @param graph The input graph.
     */
    VertexColoring(const NetworKit::Graph &graph);

    /**
     * Return the coloring found by the algorithm.
     *
     * @return a map from nodes to colors.
     */
    const std::map<NetworKit::node, int>& getColoring() const;

protected:
    const std::optional<NetworKit::Graph> graph;
    std::map<NetworKit::node, int> colors;
};

} /* namespace Koala */
