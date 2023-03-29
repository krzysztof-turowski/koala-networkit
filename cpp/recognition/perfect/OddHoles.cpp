/*
 * OddHoles.cpp
 *
 *  Created on: 05.11.2022
 *      Author: Adrian Siwiec
 *      Ported by: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include <cassert>
#include <ranges>

#include <networkit/linkprediction/NeighborhoodUtility.hpp>

#include <graph/GraphTools.hpp>
#include <traversal/BFS.hpp>
#include <traversal/DFS.hpp>
#include <traversal/PathInplace.hpp>
#include <recognition/PerfectGraphRecognition.hpp>

namespace Koala {

bool is_path(const NetworKit::Graph &graph, const std::vector<NetworKit::node> &path) {
    for (unsigned i = 0; i < path.size() - 1; i++) {
        if (!graph.hasEdge(path[i], path[i + 1])) {
            return false;
        }
        for (unsigned j = i + 2; j < path.size(); j++) {
            if (graph.hasEdge(path[i], path[j])) {
                return false;
            }
        }
    }
    return true;
}

bool PerfectGraphRecognition::containsOddHole(const NetworKit::Graph &graph) {
    std::vector<NetworKit::node> path;
    return Koala::Traversal::NextPathInplace(
        graph, std::numeric_limits<NetworKit::count>::max(), path,
        Koala::Traversal::PathInplaceMode::INDUCED_ODD_HOLE);
}

bool PerfectGraphRecognition::containsHole(
        const NetworKit::Graph &graph, NetworKit::count length) {
    if (length <= 3) {
        return false;
    }
    std::vector<NetworKit::node> path;
    return Koala::Traversal::NextPathInplace(
        graph, length, path, Koala::Traversal::PathInplaceMode::INDUCED_CYCLE);
}

bool PerfectGraphRecognition::containsT1(const NetworKit::Graph &graph) {
    return containsHole(graph, 5);
}

bool PerfectGraphRecognition::containsT2(const NetworKit::Graph &graph) {
    for (const auto &v1 : graph.nodeRange()) {
        for (const auto &v2 : graph.neighborRange(v1)) {
            for (const auto &v3 : graph.neighborRange(v2)) {
                if (v3 == v1) {
                    continue;
                }
                for (const auto &v4 : graph.neighborRange(v3)) {
                    if (v4 == v1 || v4 == v2) {
                        continue;
                    }
                    if (!is_path(graph, std::vector<NetworKit::node>{v1, v2, v3, v4})) {
                        continue;
                    }
                    auto auxiliary_components = getAuxiliaryComponents(graph, {v1, v2, v4});
                    for (auto X : auxiliary_components) {
                        if (X.empty()) {
                            continue;
                        }
                        bool path = Koala::Traversal::BFS(graph, v1, v4, [&](auto v) {
                            if (v == v2 || v == v3) {
                                return false;
                            }
                            if (graph.hasEdge(v, v2) || graph.hasEdge(v, v3)) {
                                return false;
                            }
                            if (isComplete(graph, X, v)) {
                                return false;
                            }
                            return true;
                        });
                        if (path) {
                            return true;
                        }
                    }
                }
            }
        }
    }
    return false;
}

bool is_t3(const NetworKit::Graph &graph, const std::vector<NetworKit::node> &V,
        const std::vector<NetworKit::node> &P, const std::vector<NetworKit::node> &X) {
    if (V.size() != 6 || P.empty() || X.empty()) {
        return false;
    }
    if (std::set<NetworKit::node>(V.begin(), V.end()).size() != 6) {
        return false;
    }
    std::initializer_list<std::pair<NetworKit::node, NetworKit::node>> edges{
        {0, 1}, {0, 3}, {1, 2}, {2, 3}, {2, 4}, {3, 5}};
    for (const auto &[u, v] : edges) {
        if (!graph.hasEdge(V[u], V[v])) {
            return false;
        }
    }
    std::initializer_list<std::pair<NetworKit::node, NetworKit::node>> non_edges{
        {0, 2}, {0, 4}, {0, 5}, {1, 3}, {1, 4}, {1, 5}, {3, 4}};
    for (const auto &[u, v] : non_edges) {
        if (graph.hasEdge(V[u], V[v])) {
            return false;
        }
    }

    std::vector<NetworKit::node> X_sorted(X.begin(), X.end());
    std::sort(X_sorted.begin(), X_sorted.end());
    auto auxiliary_components = PerfectGraphRecognition::getAuxiliaryComponents(
        graph, {V[0], V[1], V[4]});
    bool is_x_anticomponent = std::any_of(
        auxiliary_components.begin(), auxiliary_components.end(), [&](auto ac) {
            std::sort(ac.begin(), ac.end());
            return ac == X_sorted;
    });
    if (!is_x_anticomponent) {
        return false;
    }

    if (PerfectGraphRecognition::isComplete(graph, X, V[2])) {
        return false;
    }
    if (PerfectGraphRecognition::isComplete(graph, X, V[3])) {
        return false;
    }
    if ((P[0] != V[4] || P.back() != V[5]) && (P[0] != V[5] || P.back() != V[4])) {
        return false;
    }
    if (!is_path(graph, P)) {
        return false;
    }
    for (unsigned i = 1; i < P.size() - 1; i++) {
        if (std::ranges::any_of(std::views::iota(0, 4), [&](auto j) { return V[j] == P[i]; })) {
            return false;
        }
        if (std::any_of(X.begin(), X.end(), [&](auto x) { return x == P[i]; })) {
            return false;
        }
        if (PerfectGraphRecognition::isComplete(graph, X, P[i])) {
            return false;
        }
        if (graph.hasEdge(V[0], P[i]) || graph.hasEdge(V[1], P[i])) {
            return false;
        }
    }
    return !graph.hasEdge(V[4], V[5]) || !PerfectGraphRecognition::isComplete(graph, X, V[5]);
}

bool PerfectGraphRecognition::containsT3(const NetworKit::Graph &graph) {
    for (const auto &v1 : graph.nodeRange()) {
        for (const auto &v2 : graph.neighborRange(v1)) {
            for (const auto &v5 : graph.nodeRange()) {
                if (v5 == v1 || v5 == v2 || graph.hasEdge(v5, v1) || graph.hasEdge(v5, v2)) {
                    continue;
                }
                for (const auto &X : PerfectGraphRecognition::getAuxiliaryComponents(
                        graph, {v1, v2, v5})) {
                    if (X.empty()) {
                        continue;
                    }
                    std::set<NetworKit::node> Fprim;
                    Koala::Traversal::DFSFrom(
                        graph, v5, [&](auto v) { Fprim.insert(v); },
                        [&](auto v) {
                            if (graph.hasEdge(v1, v) || graph.hasEdge(v2, v)) {
                                return false;
                            }
                            if (isComplete(graph, X, v)) {
                                return false;
                            }
                            return true;
                    });
                    std::set<NetworKit::node> F(Fprim.begin(), Fprim.end());
                    for (const auto &fp : Fprim) {
                        for (const auto &v : graph.neighborRange(fp)) {
                            if (!F.count(v) && isComplete(graph, X, v)) {
                                if (!graph.hasEdge(v, v1) && !graph.hasEdge(v, v2)
                                        && !graph.hasEdge(v, v5)) {
                                    F.insert(v);
                                }
                            }
                        }
                    }
                    for (const auto &v4 : graph.neighborRange(v1)) {
                        if (graph.hasEdge(v4, v2) || graph.hasEdge(v4, v5)) {
                            continue;
                        }
                        auto it6 = std::find_if(F.begin(), F.end(), [&](auto f) {
                            return graph.hasEdge(v4, f);
                        });
                        if (it6 == F.end()) {
                            continue;
                        }
                        auto v6 = *it6;

                        bool v4_has_nonneighbour_in_x = std::any_of(
                            X.begin(), X.end(), [&](auto x) {
                                return !graph.hasEdge(v4, x);
                        });
                        if (!v4_has_nonneighbour_in_x) {
                            continue;
                        }

                        for (const auto &v3 : graph.neighborRange(v2)) {
                            if (!graph.hasEdge(v3, v4) || !graph.hasEdge(v3, v5)
                                    || graph.hasEdge(v3, v1)) {
                                continue;
                            }
                            bool v3_has_nonneighbour_in_x = std::any_of(
                                X.begin(), X.end(), [&](auto x) {
                                    return !graph.hasEdge(v3, x);
                            });
                            if (!v3_has_nonneighbour_in_x) {
                                continue;
                            }
                            auto P = Koala::Traversal::BFSPath(
                                graph, v6, v5, [&](int v) { return Fprim.count(v); });
                            return true;
                        }
                    }
                }
            }
        }
    }
    return false;
}

} /* namespace Koala */
