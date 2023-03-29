/*
 * MaximumFlow.hpp
 *
 *  Created on: 29.03.2023
 *      Author: Michał Stobierski
 *      Ported by: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

#include <optional>
#include <map>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace Koala {

/**
 * @ingroup flow
 * The base class for the max flow algorithms.
 *
 */
class MaximumFlow : public NetworKit::Algorithm {

public:
    /**
     * Given an input graph, set up the greedy vertex coloring procedure.
     *
     * @param graph The input graph.
     * @param s     The source vertex.
     * @param t     The sink vertex.
     */
    MaximumFlow(const NetworKit::Graph &graph, NetworKit::node s, NetworKit::node t);

    /**
     * Return the flow size found by the algorithm.
     *
     * @return a total flow value.
     */
    int getFlowSize() const;

protected:
    const std::optional<NetworKit::Graph> graph;
    NetworKit::node source, target;
    int flow_size;
};

/**
 * @ingroup flow
 * The class for the King-Rao-Tarjan maximum flow algorithm
 */
class KingRaoTarjanMaximumFlow final : public MaximumFlow {

public:
    using MaximumFlow::MaximumFlow;

    /**
     * Execute the King-Rao-Tarjan maximum flow algorithm.
     */
    void run();
};

} /* namespace Koala */
