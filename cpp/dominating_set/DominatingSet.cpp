/*
 * DominatingSet.cpp
 *
 *  Created on: 01.07.2023
 *      Author: Piotr Kubaty
 */

#include <set>
#include <vector>

#include "dominating_set/DominatingSet.hpp"

namespace Koala {

DominatingSet::DominatingSet(NetworKit::Graph &graph) : graph(std::make_optional(graph)) { }

const std::set<NetworKit::node>& DominatingSet::getDominatingSet() const {
    assureFinished();
    return dominating_set;
}

void DominatingSet::check() const {
    assureFinished();
    std::vector<bool> dominated(graph->upperNodeIdBound());
    for (const auto &u : dominating_set) {
        dominated[u] = true;
        graph->forNeighborsOf(u, [&dominated](NetworKit::node v) {
            dominated[v] = true;
        });
    }
    graph->forNodes([&dominated](NetworKit::node u) {
        assert(dominated[u]);
    });
}

}  /* namespace Koala */
