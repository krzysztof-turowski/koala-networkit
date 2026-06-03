/*
 * MaximumFlow.cpp
 *
 *  Created on: 29.03.2023
 *      Author: Michał Stobierski
 *      Ported by: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include "flow/MaximumFlow.hpp"

namespace Koala {

MaximumFlow::MaximumFlow(NetworKit::Graph &graph, NetworKit::node s, NetworKit::node t)
    : graph(graph), source(s), target(t) { }

int MaximumFlow::getFlowSize() const {
    assureFinished();
    return flow_size;
}

}  // namespace Koala
