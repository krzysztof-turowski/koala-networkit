/*
 * PathInplace.cpp
 *
 *  Created on: 22.12.2022
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include <set>

#include <traversal/PathInplace.hpp>

namespace Koala {

namespace Traversal {

bool check_last_vertex(
        const NetworKit::Graph &graph, std::vector<NetworKit::node> &path, PathInplaceMode mode) {
    switch (mode) {
        case PathInplaceMode::INDUCED_CYCLE:
            if (path.size() <= 2) {
                return false;
            }
            break;
        case PathInplaceMode::INDUCED_ODD_HOLE:
            if (path.size() <= 3 || path.size() % 2 == 0) {
                return false;
            }
            break;
        default:
            if (path.size() <= 1) {
                return false;
            }
            break;
    }
    if (std::set<NetworKit::node>(path.begin(), path.end()).size() < path.size()) {
        return false;
    }
    const auto &previous = path[path.size() - 2], &last = path.back();
    for (const auto &v : path) {
        if (v == previous) {
            if (!graph.hasEdge(previous, last)) {
                return false;
            }
        } else if (v == path[0]) {
            if (mode == PathInplaceMode::INDUCED_CYCLE || mode == PathInplaceMode::INDUCED_ODD_HOLE) {
                if (!graph.hasEdge(v, last)) {
                    return false;
                }
            } else if (mode == PathInplaceMode::INDUCED_PATH) {
                if (graph.hasEdge(v, last)) {
                    return false;
                }
            }
        } else if (v != last) {
            if (graph.hasEdge(v, last)) {
                return false;
            }
        }
    }
    if (mode == PathInplaceMode::INDUCED_ODD_HOLE) {
        for (int i = 0; i < path.size() - 1; i++) {
            for (int j = i + 2; j < path.size(); j++) {
                if (i == 0 && j == path.size() - 1) {
                    assert(graph.hasEdge(path[i], path[j]));
                } else {
                    assert(!graph.hasEdge(path[i], path[j]));
                }
            }
            assert(graph.hasEdge(path[i], path[i + 1]));
        }
    }
    return true;
}

bool NextPathInplace(
        const NetworKit::Graph &graph, NetworKit::count length, std::vector<NetworKit::node> &path,
        PathInplaceMode mode) {
    auto get_next_neighbor = [&](NetworKit::node &u, NetworKit::node &v) {
        return graph.getIthNeighbor(u, graph.indexOfNeighbor(u, v) + 1);
    };
    length = std::min(length, graph.numberOfNodes());
    path.reserve(length);
    if (path.empty()) {
        path.push_back(*graph.nodeRange().begin());
    }
    while (true) {
        if (path.back() == NetworKit::none) {
            path.pop_back();
            if (path.size() == 1) {
                auto next = ++NetworKit::Graph::NodeIterator(&graph, path[0]);
                if (next == graph.nodeRange().end()) {
                    return false;
                }
                path.back() = *next;
            } else {
                path.back() = get_next_neighbor(path[path.size() - 2], path.back());
            }
            continue;
        }
        if (path.size() < length) {
            if (path.size() > 1) {
                while (path.back() != NetworKit::none) {
                    if (mode == PathInplaceMode::INDUCED_PATH) {
                        if (check_last_vertex(graph, path, PathInplaceMode::INDUCED_PATH)) {
                            break;
                        }
                    } else if (mode == PathInplaceMode::INDUCED_CYCLE) {
                        if (path[0] < path.back() && check_last_vertex(graph, path, PathInplaceMode::INDUCED_PATH)) {
                            break;
                        }
                    } else if (mode == PathInplaceMode::INDUCED_ODD_HOLE) {
                        if (path[0] < path.back() && check_last_vertex(graph, path, PathInplaceMode::INDUCED_ODD_HOLE)) {
                            return true;
                        }
                        if (path[0] < path.back() && check_last_vertex(graph, path, PathInplaceMode::INDUCED_PATH)) {
                            break;
                        }
                    }
                    path.back() = get_next_neighbor(path[path.size() - 2], path.back());
                }
                if (path.back() == NetworKit::none) {
                    continue;
                }
            }
            path.push_back(graph.getIthNeighbor(path.back(), 0));
            if ((mode == PathInplaceMode::INDUCED_ODD_HOLE || path.size() == length) && check_last_vertex(graph, path, mode)) {
                return true;
            } else {
                continue;
            }
        }
        do {
            path.back() = get_next_neighbor(path[path.size() - 2], path.back());
        } while (path.back() != NetworKit::none && !check_last_vertex(graph, path, mode));
        if (path.back() != NetworKit::none) {
            break;
        }
    }
    return true;
}

} // namespace Traversal

} // namespace Koala
