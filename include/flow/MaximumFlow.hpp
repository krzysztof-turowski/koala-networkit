/*
 * MaximumFlow.hpp
 *
 *  Created on: 29.03.2023
 *      Author: Michał Stobierski
 *      Ported by: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/AdjListGraph.hpp>

namespace Koala {

/**
 * @ingroup flow
 * Base class for maximum-flow algorithms.
 *
 * Implementations own an internal copy of the input graph and report the value of an
 * s-t maximum flow after run() completes.
 */
class MaximumFlow : public NetworKit::Algorithm {
 public:
    /**
     * Construct a maximum-flow instance.
     *
     * @param graph The input graph.
     * @param s     The source vertex.
     * @param t     The sink vertex.
     */
    MaximumFlow(NetworKit::Graph &graph, NetworKit::node s, NetworKit::node t);

    /**
     * Return the flow size found by the algorithm.
     *
     * @return The value of the computed s-t flow.
     */
    int getFlowSize() const;

 protected:
    NetworKit::Graph graph;
    NetworKit::node source, target;
    int flow_size;
};

}  // namespace Koala
