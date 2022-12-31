/*
 * SimpleIndependentSet.hpp
 *
 *  Created on: 30.12.2022
 *      Author: Artur Salawa
 */

#pragma once

#include <map>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace Koala {

/**
 * @ingroup independentSet
 * The base class for the independent set problem algorithms.
 *
 */
class SimpleIndependentSet : public NetworKit::Algorithm {

public:
    /**
     * Given an input graph, set up the independent set problem procedure.
     *
     * @param graph The input graph.
     */
    SimpleIndependentSet(const NetworKit::Graph &graph);

    /**
     * Return the independent set found by the algorithm.
     *
     * @return a map from nodes to (true <=> belongs to the independent set).
     */
    const std::map<NetworKit::node, bool>& getIndependentSet() const;

protected:
    const std::optional<NetworKit::Graph> graph;
    std::map<NetworKit::node, bool> independentSet;
};

/**
 * @ingroup independentSet
 * The class for the brute force algorithm.
 */
class BruteForceIndependentSet final : public SimpleIndependentSet {

public:
    using SimpleIndependentSet::SimpleIndependentSet;

    /**
     * Execute the brute force algorithm.
     */
    void run();
};

} /* namespace Koala */
