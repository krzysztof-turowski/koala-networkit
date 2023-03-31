/*
 * VertexColoring.cpp
 *
 *  Created on: 30.03.2023
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include <coloring/VertexColoring.hpp>

namespace Koala {

VertexColoring::VertexColoring(
        const NetworKit::Graph &graph) : graph(std::make_optional(graph)) { }

const std::map<NetworKit::node, int>& VertexColoring::getColoring() const {
    assureFinished();
    return colors;
}

} /* namespace Koala */
