/*
 * BFS.hpp
 *
 *  Created on: 12.11.2021
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

#include <map>
#include <queue>
#include <vector>

#include <networkit/graph/Graph.hpp>

namespace Koala {

namespace Traversal {

template <typename Predicate>
bool BFS(
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

template <typename Predicate>
std::vector<NetworKit::node> BFSPath(
        const NetworKit::Graph &G, NetworKit::node source, NetworKit::node target,
        Predicate predicate) {
    std::queue<NetworKit::node> Q({source});
    std::map<NetworKit::node, NetworKit::node> parent;
    parent[source] = source;
    while (!Q.empty()) {
        const auto u = Q.front();
        Q.pop();
        if (u == target) {
            break;
        }
        G.forNeighborsOf(u, [&](NetworKit::node v) {
            if (!parent.count(v) && (predicate(v) || v == target)) {
                Q.push(v);
                parent[v] = u;
            }
        });
    }
    if (!parent.count(target)) {
        return std::vector<NetworKit::node>();
    }

    auto v = target;
    std::vector<NetworKit::node> out({target});
    while (parent[v] != v) {
        v = parent[v];
        out.push_back(v);
    }
    std::reverse(out.begin(), out.end());
    return out;
}

} // namespace Traversal

} // namespace Koala
