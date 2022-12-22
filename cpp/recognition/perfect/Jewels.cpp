/*
 * Jewels.cpp
 *
 *  Created on: 04.11.2022
 *      Author: Adrian Siwiec
 *      Ported by: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include <cassert>
#include <set>
#include <vector>

#include <networkit/linkprediction/NeighborhoodUtility.hpp>

#include <traversal/BFS.hpp>
#include <traversal/PathInplace.hpp>
#include <recognition/PerfectGraphRecognition.hpp>

#include "../other/jewels.h"

namespace Koala {

bool is_jewel(const NetworKit::Graph &graph, const std::vector<NetworKit::node> &vertices) {
    if (vertices.size() != 5) {
        return false;
    }
    std::set<NetworKit::node> vertices_set(vertices.begin(), vertices.end());
    if (vertices_set.size() != 5) {
        return false;
    }
    std::initializer_list<std::pair<NetworKit::node, NetworKit::node>> non_edges = {
        {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}};
    for (auto &[u, v] : non_edges) {
        if (!graph.hasEdge(vertices[u], vertices[v])) {
            return false;
        }
    }
    std::initializer_list<std::pair<NetworKit::node, NetworKit::node>> edges = {
        {0, 2}, {1, 3}, {0, 3}};
    for (auto &[u, v] : edges) {
        if (graph.hasEdge(vertices[u], vertices[v])) {
            return false;
        }
    }

    std::set<NetworKit::node> vertices124{vertices[1], vertices[2], vertices[4]};
    auto non_neighbours = [&](NetworKit::node v) {
        if (vertices124.count(v)) {
            return false;
        }
        for (const auto &u : graph.neighborRange(v)) {
            if (vertices124.count(u)) {
                return false;
            }
        }
        return true;
    };
    return Koala::Traversal::BFS(graph, vertices[0], vertices[3], non_neighbours);
}

// TODO(kturowski): temporary check
void check_jewel(const NetworKit::Graph &graph, Graph &G, std::vector<NetworKit::node> &path) {
    vec<int> v(path.begin(), path.end());
    assert(is_jewel(graph, path) == isJewel(G, v));
}

bool PerfectGraphRecognition::containsJewel(const NetworKit::Graph &graph) {
    std::vector<NetworKit::node> path;
    Graph G(graph);
    while (Koala::Traversal::NextPathInplace(
            graph, 4, path, Koala::Traversal::PathInplaceMode::INDUCED_PATH)) {
        for (auto v5 : NetworKit::NeighborhoodUtility::getCommonNeighbors(
                graph, path[0], path.back())) {
            path.push_back(v5);
            check_jewel(graph, G, path);
            if (is_jewel(graph, path)) {
                return true;
            }
            path.pop_back();
        }
    }
    assert(!::containsJewelNaive(G));
    return false;
}

} /* namespace Koala */
