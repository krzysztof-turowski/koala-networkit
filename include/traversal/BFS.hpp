/*
 * BFS.hpp
 *
 *  Created on: 12.11.2021
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

#include <queue>
#include <vector>

#include <networkit/graph/Graph.hpp>

namespace Koala {

namespace Traversal {

template <typename Predicate>
bool BFSwithPredicate(
        const NetworKit::Graph &G, NetworKit::node source, NetworKit::node target,
        Predicate predicate) {
    std::queue<NetworKit::node> Q({source});
    std::vector<bool> marked(G.upperNodeIdBound());
    marked[source] = true;
    while (!Q.empty()) {
        const auto u = Q.front();
        Q.pop();
        if (u == target) {
            return true;
        }
        G.forNeighborsOf(u, [&](NetworKit::node v) {
            if (!marked[v] && (predicate(v) || v == target)) {
                Q.push(v);
                marked[v] = true;
            }
        });
    }
    return false;
}

} // namespace Traversal

} // namespace Koala
