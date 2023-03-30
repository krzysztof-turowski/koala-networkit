/*
 * Pyramids.cpp
 *
 *  Created on: 20.12.2022
 *      Author: Adrian Siwiec
 *      Ported by: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include <algorithm>
#include <cassert>
#include <ranges>
#include <set>
#include <vector>

#include <traversal/BFS.hpp>
#include <recognition/PerfectGraphRecognition.hpp>

namespace Koala {

void nextTupleInPlace(std::vector<NetworKit::node> &vertices, NetworKit::count max) {
    vertices[0]++;
    for (unsigned i = 0; i < vertices.size() && vertices[i] >= max; i++) {
        vertices[i] = 0;
        if (i + 1 < vertices.size()) {
            vertices[i + 1]++;
        }
    }
}

std::vector<NetworKit::node> nextTuple(
        std::vector<NetworKit::node> &vertices, NetworKit::count max) {
    nextTupleInPlace(vertices, max);
    return vertices;
}

std::vector<std::vector<NetworKit::node>> generateTuples(
        NetworKit::count size, NetworKit::count max) {
    std::vector<std::vector<NetworKit::node>> out;
    auto current = std::vector<NetworKit::node>(size);
    do {
        out.push_back(current);
        current = nextTuple(current, max);
    } while (std::any_of(current.cbegin(), current.cend(), [](auto i){ return i != 0; }));
    return out;
}

bool check_prerequisites(
        const NetworKit::Graph &graph, NetworKit::node a, const std::vector<NetworKit::node> &b,
        const std::vector<NetworKit::node> &s) {
    // We assume b is a Triangle and s is an EmptyStarTriangle
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            if (i == j) {
                continue;
            }
            if (b[i] == s[j]) {
                return false;
            }
            if ((b[j] != s[j] && graph.hasEdge(b[i], s[j])) || graph.hasEdge(s[i], s[j])) {
                return false;
            }
        }
    }

    bool aAdjB = false;
    for (int i = 0; i < 3; i++) {
        if (!graph.hasEdge(a, s[i])) {
            return false;
        }
        if (graph.hasEdge(a, b[i])) {
            if (aAdjB || b[i] != s[i]) {
                return false;
            }
            aAdjB = true;
        }
    }
    return true;
}

std::vector<std::vector<NetworKit::node>> get_all_triangles(const NetworKit::Graph &graph) {
    std::vector<std::vector<NetworKit::node>> out;
    for (auto i : graph.nodeRange()) {
        for (auto j : graph.neighborRange(i)) {
            for (auto k : graph.neighborRange(j)) {
                if (i < j && j < k && graph.hasEdge(i, k)) {
                    out.emplace_back(std::vector<NetworKit::node>{i, j, k});
                }
            }
        }
    }
    return out;
}

auto get_all_empty_star_triangles(const NetworKit::Graph &graph) {
    std::vector<std::pair<NetworKit::node, std::vector<NetworKit::node>>> out;
    for (const auto &a : graph.nodeRange()) {
        for (const auto &s1 : graph.neighborRange(a)) {
            for (const auto &s2 : graph.neighborRange(a)) {
                if (s2 == s1 || graph.hasEdge(s1, s2)) {
                    continue;
                }
                for (const auto &s3 : graph.neighborRange(a)) {
                    if (s3 == s2 || s3 == s1 || graph.hasEdge(s1, s3) || graph.hasEdge(s2, s3)) {
                        continue;
                    }
                    out.push_back(std::make_pair(a, std::vector<NetworKit::node>{s1, s2, s3}));
                }
            }
        }
    }
    return out;
}

bool vectors_cut_empty(auto a_begin, auto a_end, auto b_begin, auto b_end) {
    std::set<int> a(a_begin, a_end);
    for (auto it = b_begin; it != b_end; ++it) {
        if (a.count(*it) > 0) {
            return false;
        }
    }
    return true;
}

bool no_edges_between_vectors(
        const NetworKit::Graph &graph, auto a_begin, auto a_end, auto b_begin, auto b_end) {
    std::set<int> a_neighborhood;
    for (auto it = a_begin; it != a_end; ++it) {
        for (int i : graph.neighborRange(*it)) {
            a_neighborhood.insert(i);
        }
    }
    for (auto it = b_begin; it != b_end; ++it) {
        if (a_neighborhood.count(*it) > 0) {
          return false;
        }
    }
    return true;
}

auto calculate_p_paths(
        const NetworKit::Graph &graph, const std::vector<NetworKit::node> &b,
        const std::vector<NetworKit::node> &s, const std::set<NetworKit::node> &M) {
    // S[a][b][c] = S_a(b)[c], c-th vertex of the Sa(b) path
    std::vector<std::vector<NetworKit::node>> S[3], T[3];
    for (int i = 0; i < 3; i++) {
        S[i].resize(graph.upperNodeIdBound()), T[i].resize(graph.upperNodeIdBound());
        for (const auto &m : graph.nodeRange()) {
            auto non_neighbours = [&](NetworKit::node v) {
                if (v == s[i] || v == m) {
                    return true;
                }
                if (M.count(v)) {
                    return false;
                }
                for (int j = 0; j < 3; j++) {
                    if ((i != j) && (graph.hasEdge(s[j], v) || graph.hasEdge(v, b[j]))) {
                        return false;
                    }
                }
                return true;
            };
            S[i][m] = Koala::Traversal::BFSPath(graph, s[i], m, non_neighbours);
            T[i][m] = Koala::Traversal::BFSPath(graph, m, b[i], non_neighbours);
        }
    }

    std::vector<std::vector<std::vector<NetworKit::node>>> P(3);
    for (int i = 0; i < 3; i++) {
        P[i].resize(graph.upperNodeIdBound());
        if (s[i] == b[i]) {
            P[i][b[i]] = std::vector<NetworKit::node>{b[i]};
        } else {
            for (const auto &m : graph.nodeRange()) {
                if ((S[i][m].empty() && s[i] != m) || (T[i][m].empty() && m != b[i])) {
                    continue;
                }
                bool noPath = !M.count(m) && std::ranges::any_of(
                    std::views::iota(0, 3), [&](auto j) {
                        return (i != j) && (graph.hasEdge(s[j], m) || graph.hasEdge(m, b[j]));
                });
                if (noPath) {
                    continue;
                }
                auto &sim = S[i][m];
                auto &tim = T[i][m];
                if (!vectors_cut_empty(sim.begin(), sim.end() - 1, tim.begin() + 1, tim.end())) {
                    continue;
                }
                if (!no_edges_between_vectors(
                        graph, sim.begin(), sim.end() - 1, tim.begin() + 1, tim.end())) {
                    continue;
                }
                P[i][m].insert(P[i][m].end(), sim.begin(), sim.end());
                P[i][m].insert(P[i][m].end(), tim.begin() + 1, tim.end());
            }
        }
    }
    return P;
}

bool is_pyramid(
        const NetworKit::Graph &graph, NetworKit::node a, const std::vector<NetworKit::node> &b,
        const std::vector<std::vector<NetworKit::node>> &P) {
    if (b.size() != 3 || P.size() != 3) {
        return false;
    }
    for (const auto &p : P) {
        if (p.empty()) {
            return false;
        }
    }
    for (int i = 0; i < 3; i++) {
        for (int j = i + 1; j < 3; j++) {
            if (b[i] == b[j] || !graph.hasEdge(b[i], b[j])) {
                return false;
            }
            if (P[i][0] == P[j][0] || graph.hasEdge(P[i][0], P[j][0])) {
                return false;
            }
        }
    }
    for (int i = 0; i < 3; i++) {
        if (!graph.hasEdge(a, P[i][0]) || P[i].back() != b[i]) {
            return false;
        }
    }
    for (const auto &p : P) {
        for (unsigned i = 1; i < p.size(); i++) {
            if (!graph.hasEdge(p[i - 1], p[i])) {
                return false;
            }
        }
    }
    for (int i = 0; i < 3; i++) {
        for (int j = i + 1; j < 3; j++) {
            int edges = 0;
            for (int a : P[i]) {
                for (int b : P[j]) {
                    if (graph.hasEdge(a, b)) {
                        edges++;
                    }
                }
            }
            if (edges != 1) {
                return false;
            }
        }
    }
    return std::count_if(b.cbegin(), b.cend(), [&](auto v) { return graph.hasEdge(a, v); }) <= 1;
}

bool PerfectGraphRecognition::containsPyramid(const NetworKit::Graph &graph) {
    auto triangles = get_all_triangles(graph);
    auto emptyStarTriangles = get_all_empty_star_triangles(graph);
    for (const auto &b : triangles) {
        for (const auto &[a, s] : emptyStarTriangles) {
            if (!check_prerequisites(graph, a, b, s)) {
                continue;
            }
            std::set<NetworKit::node> M(s.begin(), s.end());
            M.insert(b.begin(), b.end());
            auto P = calculate_p_paths(graph, b, s, M);

            std::vector<std::pair<int, int>> pairs{{0, 1}, {1, 2}, {2, 0}};
            std::set<std::pair<NetworKit::node, NetworKit::node>> goodPairs[3];
            for (const auto &[u, v] : pairs) {
                for (const auto &m1 : graph.nodeRange()) {
                    if (P[u][m1].empty()) {
                        continue;
                    }

                    std::set<NetworKit::node> color;
                    for (auto i : P[u][m1]) {
                        if (M.count(i)) {
                            continue;
                        }
                        color.insert(i);
                        for (const auto &j : graph.neighborRange(i)) {
                            color.insert(j);
                        }
                    }

                    for (const auto &m2 : graph.nodeRange()) {
                        if (P[v][m2].empty()) {
                            continue;
                        }
                        bool found = std::none_of(P[v][m2].cbegin(), P[v][m2].cend(), [&](auto i) {
                            return color.count(i);
                        });
                        if (found) {
                            goodPairs[u].emplace(std::make_pair(m1, m2));
                        }
                    }
                }
            }  // (i, j) good pairs completed

            auto triples = generateTuples(3, graph.numberOfNodes());
            for (const auto &triple : triples) {
                bool found = std::all_of(pairs.cbegin(), pairs.cend(), [&](const auto &pair) {
                    const auto& [u, v] = pair;
                    return goodPairs[u].count(std::make_pair(triple[u], triple[v]));
                });
                if (found) {
                  std::vector<std::vector<NetworKit::node>> paths = {
                      P[0][triple[0]], P[1][triple[1]], P[2][triple[2]]
                  };
                  return true;
                }
            }
        }
    }
    return false;
}

} /* namespace Koala */
