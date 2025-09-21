/*
 * GraphTools.hpp
 *
 *  Created on: 12.11.2021
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

#include <networkit/graph/Graph.hpp>

namespace Koala {

namespace GraphTools {

NetworKit::Graph toComplement(const NetworKit::Graph&);

NetworKit::Graph convertDirectedGraphToUndirected(NetworKit::Graph&, bool weighted = false);
NetworKit::Graph convertUndirectedGraphToDirected(NetworKit::Graph&, bool weighted = false);

void printGraph(const NetworKit::Graph&);

}  // namespace GraphTools

}  // namespace Koala
