/*
 * DominatingSet.hpp
 *
 *  Created on: 01.07.2023
 *      Author: Piotr Kubaty
 */

#pragma once

#include <optional>
#include <set>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace Koala {

/**
 * @ingroup dominating_set
 * The base class for the dominating set algorithms.
 *
 */
class DominatingSet : public NetworKit::Algorithm {
 public:
    /**
     * Given an input graph, set up the vertex coloring procedure.
     *
     * @param graph The input graph.
     */
    explicit DominatingSet(NetworKit::Graph &graph);

    /**
     * Return the dominating set found by the algorithm.
     *
     * @return a set with vertices in the dominating set.
     */
    const std::set<NetworKit::node>& getDominatingSet() const;

    /**
     * Verify the result found by the algorithm.
     */
    void check() const;

 protected:
    std::optional<NetworKit::Graph> graph;
    std::set<NetworKit::node> dominating_set;
};

}  /* namespace Koala */
