/*
 * IndependentSet.hpp
 *
 *  Created on: 30.12.2022
 *      Author: Artur Salawa
 */

#pragma once

#include <map>
#include <optional>
#include <set>
#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace Koala {

/**
 * @ingroup independent_set
 * The base class for independent set algorithms.
 *
 */
class IndependentSet : public NetworKit::Algorithm {
 public:
    /**
     * Given an input graph, set up the independent set procedure.
     *
     * @param graph The input graph.
     */
    explicit IndependentSet(const NetworKit::Graph &graph);

    /**
     * Return the independent set found by the algorithm.
     *
     * @return a set of nodes that form an independent set.
     */
    const std::set<NetworKit::node>& getIndependentSet() const;

    /**
     * Execute the maximum independent set finding procedure.
     */
    virtual void run() = 0;

    /**
     * Verify the result found by the algorithm.
     */
    void check() const;

 protected:
    std::optional<NetworKit::Graph> graph;
    std::set<NetworKit::node> independentSet;
};

} /* namespace Koala */
