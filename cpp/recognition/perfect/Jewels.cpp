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

namespace Koala {

bool is_jewel(const NetworKit::Graph &graph, const std::vector<NetworKit::node> &V) {
    if (V.size() != 5) {
        return false;
    }
    std::set<NetworKit::node> vertices_set(V.begin(), V.end());
    if (vertices_set.size() != 5) {
        return false;
    }
    std::initializer_list<std::pair<NetworKit::node, NetworKit::node>> non_edges{
        {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}};
    for (const auto &[u, v] : non_edges) {
        if (!graph.hasEdge(V[u], V[v])) {
            return false;
        }
    }
    std::initializer_list<std::pair<NetworKit::node, NetworKit::node>> edges{
        {0, 2}, {1, 3}, {0, 3}};
    for (const auto &[u, v] : edges) {
        if (graph.hasEdge(V[u], V[v])) {
            return false;
        }
    }

    std::set<NetworKit::node> vertices124{V[1], V[2], V[4]};
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
    return Koala::Traversal::BFS(graph, V[0], V[3], non_neighbours);
}

bool PerfectGraphRecognition::contains_jewel(const NetworKit::Graph &graph) {
    std::vector<NetworKit::node> P;
    while (Koala::Traversal::NextPathInplace(
            graph, 4, P, Koala::Traversal::PathInplaceMode::INDUCED_PATH)) {
        for (const auto &v5 : NetworKit::NeighborhoodUtility::getCommonNeighbors(
                graph, P[0], P.back())) {
            P.push_back(v5);
            if (is_jewel(graph, P)) {
                return true;
            }
            P.pop_back();
        }
    }
    return false;
}

} /* namespace Koala */
