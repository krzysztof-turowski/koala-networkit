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
    std::vector<bool> marked(G.upperNodeIdBound());
    std::queue<NetworKit::node> Q({source});
    while (!Q.empty()) {
        const auto u = Q.front();
        Q.pop();
        if (u == target) {
            return true;
        }
        if (predicate(u)) {
            G.forNeighborsOf(u, [&](NetworKit::node v) {
                if (!marked[v]) {
                    Q.push(v);
                    marked[v] = true;
                }
            });
        }
    }
    return false;
}

} // namespace Traversal

} // namespace Koala
