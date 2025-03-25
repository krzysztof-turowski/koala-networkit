/*
 * MaxCut.cpp
 *
 * Solution for the Max-Cut.
 * Created on: 13.05.2024
 * Author: Michał Miziołek
 */

#include <ranges>

#include <max_cut/MaxCut.hpp>

namespace Koala {

MaxCut::MaxCut(NetworKit::Graph &graphInput)
    : graph(std::make_optional(graphInput)) {}

int MaxCut::getMaxCutValue() const {
    return maxCutValue;
}

const std::vector<bool>& MaxCut::getMaxCutSet() const {
    return maxCutSet;
}

double MaxCut::calculateCutValue(const std::vector<bool>& set) {
    double cutValue = 0;
    graph->forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) {
        if (set[u] != set[v]) {
            cutValue += w;
        }
    });
    return cutValue;
}

}  // namespace Koala
