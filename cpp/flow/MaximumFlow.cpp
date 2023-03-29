/*
 * MaximumFlow.cpp
 *
 *  Created on: 29.03.2023
 *      Author: Michał Stobierski
 *      Ported by: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include <iostream>

#include <flow/MaximumFlow.hpp>

#include "krt/king_rao_tarjan.hpp"

namespace Koala {

MaximumFlow::MaximumFlow(const NetworKit::Graph &graph, NetworKit::node s, NetworKit::node t)
    : graph(std::make_optional(graph)), source(s), target(t) { }

int MaximumFlow::getFlowSize() const {
    assureFinished();
    return flow_size;
}

void KingRaoTarjanMaximumFlow::run() {
    std::cout << "KingRaoTarjanMaximumFlow::run() started" << std::endl;
    flow_size = king_rao_tarjan(KRTGraph::read_from_graph(*graph, source, target));
    std::cout << "KingRaoTarjanMaximumFlow::run() finished" << std::endl;
    hasRun = true;
}

} /* namespace Koala */
