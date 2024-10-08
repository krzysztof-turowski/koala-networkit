/*
 * MinCut.cpp
 *
 * Solution for the Min-Cut.
 * Created on: 17.05.2024
 * Author: Michał Miziołek
 */

#include <ranges>

#include <min_cut/MinCut.hpp>

namespace Koala {

MinCut::MinCut(NetworKit::Graph &graphInput)
    : graph(std::make_optional(graphInput)) {}

int MinCut::getMinCutValue() const {
    return minCutValue;
}

const std::vector<bool>& MinCut::getMinCutSet() const {
    return minCutSet;
}

double MinCut::calculateCutValue(const std::vector<bool>& set) {
    double cutValue = 0;
    graph->forEdges([&](NetworKit::node u,
            NetworKit::node v, NetworKit::edgeweight w) {
        if (set[u] != set[v]) {
            cutValue += w;
        }
    });
    return cutValue;
}

}  // namespace Koala
