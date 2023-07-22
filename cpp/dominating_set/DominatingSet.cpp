/*
 * DominatingSet.cpp
 *
 *  Created on: 01.07.2023
 *      Author: Piotr Kubaty
 */

#include <dominating_set/DominatingSet.hpp>

namespace Koala {

DominatingSet::DominatingSet(NetworKit::Graph &graph) : graph(std::make_optional(graph)) { }

const std::set<NetworKit::node>& DominatingSet::getDominatingSet() const {
    assureFinished();
    return dominating_set;
}

void DominatingSet::check() const {
    assureFinished();
    std::vector<bool> dominated(graph->numberOfNodes());
    for (const auto &u : dominating_set) {
        dominated[u] = true;
        graph->forNeighborsOf(u, [&dominated](NetworKit::node v) {
            dominated[v] = true;
        });
    }
    assert(std::all_of(dominated.begin(), dominated.end(), [](bool d) { return d; }));
}

}  /* namespace Koala */
