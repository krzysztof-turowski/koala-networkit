/*
 * DFS.hpp
 *
 *  Created on: 22.12.2022
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

#include <stack>
#include <vector>

#include <networkit/graph/Graph.hpp>

namespace Koala {

namespace Traversal {

template <typename Action, typename Predicate>
void DFSFrom(
        const NetworKit::Graph &G, NetworKit::node source,
        Action action, Predicate predicate) {
    std::stack<NetworKit::node> Q({source});
    std::vector<bool> marked(G.upperNodeIdBound());
    marked[source] = true;
    while (!Q.empty()) {
        const auto u = Q.top();
        Q.pop();
        action(u);
        G.forNeighborsOf(u, [&](NetworKit::node v) {
            if (!marked[v] && predicate(v)) {
                Q.push(v);
                marked[v] = true;
            }
        });
    }
}

} // namespace Traversal

} // namespace Koala
