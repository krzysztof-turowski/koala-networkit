/*
 * Jewels.cpp
 *
 *  Created on: 04.11.2022
 *      Author: Adrian Siwiec ()
 *      Ported by: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include <cassert>
#include <set>
#include <vector>

#include <networkit/linkprediction/NeighborhoodUtility.hpp>

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
    return Koala::Traversal::BFSwithPredicate(graph, vertices[0], vertices[3], non_neighbours);
}

bool next_path_inplace(
        const NetworKit::Graph &graph, unsigned length, std::vector<NetworKit::node> &path) {
    auto get_next_neighbor = [&](NetworKit::node &u, NetworKit::node &v) {
        return graph.getIthNeighbor(u, graph.indexOfNeighbor(u, v) + 1);
    };
    auto check_last_on_path = [&](std::vector<NetworKit::node> &path) {
        if (path.size() <= 1) {
            return true;
        }
        const auto &previous = path[path.size() - 2], &last = path.back();
        if (!graph.hasEdge(previous, last)) {
            return false;
        }
        for (const auto &v : path) {
            if (v != last && v != previous && graph.hasEdge(v, last)) {
                return false;
            }
        }
        return true;
    };

    path.reserve(length);
    if (path.empty()) {
        path.push_back(*graph.nodeRange().begin());
    }
    do {
        if (path.back() == NetworKit::none) {
            path.pop_back();
            if (path.size() == 1) {
                auto next = ++NetworKit::Graph::NodeIterator(&graph, path[0]);
                if (next == graph.nodeRange().end()) {
                    return false;
                }
                path[0] = *next;
                continue;
            } else {
                path.back() = get_next_neighbor(path[path.size() - 2], path.back());
                continue;
            }
        }

        if (path.size() < length) {
            if (path.size() > 1) {
                while (path.back() != NetworKit::none && !check_last_on_path(path)) {
                    path.back() = get_next_neighbor(path[path.size() - 2], path.back());
                }
            }
            if (path.back() == NetworKit::none) {
                continue;
            }
            path.push_back(graph.getIthNeighbor(path.back(), 0));
            if (path.size() == length && check_last_on_path(path)) {
                return true;
            } else {
                continue;
            }
        }

        do {
            path.back() = get_next_neighbor(path[path.size() - 2], path.back());
        } while (path.back() != NetworKit::none && !check_last_on_path(path));
    } while (path.back() == NetworKit::none);
    
    return true;
}

// TEMPORARY
void check(const NetworKit::Graph &graph, Graph &G, std::vector<NetworKit::node> &path) {
    vec<int> v;
    for (auto& u : path) {
      v.push_back(u);
    }
    assert(is_jewel(graph, path) == isJewel(G, v));
}

bool PerfectGraphRecognition::containsJewel(const NetworKit::Graph &graph) {
    std::vector<NetworKit::node> path;
    Graph G(graph);
    while (next_path_inplace(graph, 4, path)) {
      for (auto v5 : NetworKit::NeighborhoodUtility::getCommonNeighbors(graph, path[0], path.back())) {
          path.push_back(v5);
          check(graph, G, path);
          if (is_jewel(graph, path)) {
              return true;
          }
          path.pop_back();
      }
    }
    return false;
}

} /* namespace Koala */
