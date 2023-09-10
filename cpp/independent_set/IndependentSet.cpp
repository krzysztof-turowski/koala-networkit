/*
 * IndependentSet.cpp
 *
 *  Created on: 30.12.2022
 *      Author: Artur Salawa
 */

#include <independent_set/IndependentSet.hpp>

namespace Koala {

IndependentSet::IndependentSet(const NetworKit::Graph &graph)
    : graph(std::make_optional(graph)) { }

const std::set<NetworKit::node>& IndependentSet::getIndependentSet() const {
    assureFinished();
    return independentSet;
}

void IndependentSet::check() const {
    assureFinished();
    graph->forEdges([&](NetworKit::node u, NetworKit::node v) {
        assert(!(independentSet.contains(u) && independentSet.contains(v)));
    });
}

}  /* namespace Koala */
