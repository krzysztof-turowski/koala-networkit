/*
 * IndependentSet.hpp
 *
 *  Created on: 06.04.2023
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

#include <optional>
#include <set>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace Koala {

/**
 * @ingroup independent_set
 * The base class for the independent set algorithms.
 *
 */
class IndependentSet : public NetworKit::Algorithm {
 public:
    /**
     * Given an input graph, set up the independent set procedure.
     *
     * @param graph The input graph.
     */
    IndependentSet(NetworKit::Graph &graph);

    /**
     * Return the independent set found by the algorithm.
     *
     * @return a set of nodes.
     */
    const std::set<NetworKit::node>& getIndependentSet() const;

 protected:
    std::optional<NetworKit::Graph> graph;
    std::set<NetworKit::node> out_is;
};

}  /* namespace Koala */
