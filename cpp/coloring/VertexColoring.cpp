/*
 * VertexColoring.cpp
 *
 *  Created on: 30.03.2023
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include <algorithm>
#include <cassert>
#include <map>
#include <numeric>

#include <coloring/VertexColoring.hpp>

namespace Koala {

VertexColoring::VertexColoring(const NetworKit::Graph &graph) : graph(std::make_optional(graph)) { }

const std::map<NetworKit::node, NetworKit::count>& VertexColoring::getColoring() const {
    assureFinished();
    return colors;
}

NetworKit::count VertexColoring::getMaximumColor() const {
    assureFinished();
    return std::accumulate(
        colors.begin(), colors.end(),
        static_cast<NetworKit::count>(0),
        [](NetworKit::count max_color, const auto &v) {
            return std::max(max_color, v.second);
    });
}

void VertexColoring::check() const {
    assureFinished();
    graph->forEdges([&](NetworKit::node u, NetworKit::node v) {
        assert(colors.at(u) != colors.at(v));
    });
}

} /* namespace Koala */
