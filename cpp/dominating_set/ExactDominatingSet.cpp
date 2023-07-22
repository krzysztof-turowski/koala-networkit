/*
 * ExactDominatingSet.hpp
 *
 *  Created on: 01.07.2023
 *      Author: Piotr Kubaty
 */

#include <cassert>
#include <map>
#include <ranges>
#include <tuple>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/max_cardinality_matching.hpp>

#include <dominating_set/ExactDominatingSet.hpp>

namespace Koala {

bool is_optional_dominating_set(
        const NetworKit::Graph &G, const std::set<NetworKit::node> &solution,
        std::set<NetworKit::node> bound) {
    for (const auto &u : solution) {
        bound.erase(u);
        G.forNeighborsOf(u, [&bound](NetworKit::node v) { bound.erase(v); });
    }
    return bound.empty();
}

std::vector<NetworKit::node> merge(
        const std::set<NetworKit::node> &A, const std::set<NetworKit::node> &B) {
    std::vector<NetworKit::node> merged;
    std::set_union(A.begin(), A.end(), B.begin(), B.end(), std::inserter(merged, merged.begin()));
    return merged;
}

bool ExactDominatingSet::find_small_MODS_recursive(
        const NetworKit::Graph &G, const std::vector<NetworKit::node> &V,
        NetworKit::index index, NetworKit::count size, std::set<NetworKit::node> &S) {
    if (index == V.size()) {
        return is_optional_dominating_set(G, S, bound);
    }
    if (S.size() < size) {
        S.insert(V.at(index));
        if (find_small_MODS_recursive(G, V, index + 1, size, S)) {
            return true;
        }
        S.erase(V.at(index));
    }
    if (V.size() - index <= size - S.size()) {
        return false;
    }
    return find_small_MODS_recursive(G, V, index + 1, size, S);
}

void FominKratschWoegingerDominatingSet::run() {
    hasRun = true;
    graph->forNodes([this](NetworKit::node u) {
        bound.insert(u);
        if (graph->degree(u) == 1) {
            degree[0].insert(u);
        } else if (graph->degree(u) == 2) {
            degree[1].insert(u);
        }
    });
    NetworKit::Graph G_directed(graph->numberOfNodes(), false, true);
    for (const auto &[u, v] : graph->edgeRange()) {
        G_directed.addEdge(u, v), G_directed.addEdge(v, u);
    }
    dominating_set = find_big_MODS_recursive(G_directed);
}

std::set<NetworKit::node> FominKratschWoegingerDominatingSet::find_big_MODS_recursive(
        NetworKit::Graph &G) {
    std::set<NetworKit::node> solution;
    if (!degree[0].empty()) {
        NetworKit::node u = *degree[0].begin(), unique = G.getIthNeighbor(u, 0);
        bool is_free = forget_vertex(G, u, false);
        if (is_free) {
            solution = find_big_MODS_recursive(G);
        } else {
            auto moved_unique = move_to_solution(G, unique);
            bool is_free_unique = forget_vertex(G, unique, true);
            solution = find_big_MODS_recursive(G);
            retrieve_vertex(G, unique, is_free_unique, true);
            remove_from_solution(unique, moved_unique);
        }
        retrieve_vertex(G, u, is_free, false);
        return solution;
    }
    if (!degree[1].empty()) {
        NetworKit::node v = *degree[1].begin();
        NetworKit::node u1 = G.getIthNeighbor(v, 0), u2 = G.getIthNeighbor(v, 1);
        bool is_free_v, is_free_u1, is_free_u2;
        // Branch 1
        is_free_v = forget_vertex(G, v, false);
        auto moved_1 = move_to_solution(G, u1);
        is_free_u1 = forget_vertex(G, u1, true);
        solution = find_big_MODS_recursive(G);
        required.erase(u1), retrieve_vertex(G, u1, is_free_u1, true);
        remove_from_solution(u1, moved_1);
        retrieve_vertex(G, v, is_free_v, false);
        // Branch 2
        is_free_u1 = forget_vertex(G, u1, false), is_free_u2 = forget_vertex(G, u2, false);
        auto moved_2 = move_to_solution(G, v);
        is_free_v = forget_vertex(G, v, true);
        auto other_solution = find_big_MODS_recursive(G);
        retrieve_vertex(G, v, is_free_v, true);
        remove_from_solution(v, moved_2);
        retrieve_vertex(G, u2, is_free_u2, false), retrieve_vertex(G, u1, is_free_u1, false);
        if (other_solution.size() < solution.size()) {
            solution.swap(other_solution);
        }
        // Branch 3
        is_free_v = forget_vertex(G, v, false);
        if (is_free_v) {
            other_solution = find_big_MODS_recursive(G);
        } else {
            auto moved_3 = move_to_solution(G, u2);
            is_free_u2 = forget_vertex(G, u2, true);
            other_solution = find_big_MODS_recursive(G);
            retrieve_vertex(G, u2, is_free_u2, true);
            remove_from_solution(u2, moved_3);
        }
        retrieve_vertex(G, v, is_free_v, false);
        return solution.size() < other_solution.size() ? solution : other_solution;
    }
    std::set<NetworKit::node> bound_isolated_vertices;
    for (const auto &u : G.nodeRange()) {
        if (G.degree(u) == 0 && bound.contains(u)) {
            bound.erase(u), bound_isolated_vertices.insert(u);
        }
    }
    if (!bound.empty()) {
        solution = find_MODS_for_minimum_degree_3(G);
    }
    solution.insert(required.begin(), required.end());
    solution.insert(bound_isolated_vertices.begin(), bound_isolated_vertices.end());
    bound.insert(bound_isolated_vertices.begin(), bound_isolated_vertices.end());
    return solution;
}

std::set<NetworKit::node> FominKratschWoegingerDominatingSet::find_MODS_for_minimum_degree_3(
        NetworKit::Graph &G) {
    NetworKit::Graph G_prim(G.numberOfNodes());
    for (const auto &[u, v] : G.edgeRange()) {
        if (!required.contains(u) && !required.contains(v) && !G_prim.hasEdge(u, v)) {
            G_prim.addEdge(u, v);
        }
    }
    auto possibilities = merge(free, bound);
    for (NetworKit::count i = 1; 8 * i <= 3 * possibilities.size(); i++) {
        std::set<NetworKit::node> solution;
        if (find_small_MODS_recursive(G_prim, possibilities, 0, i, solution)) {
            return solution;
        }
    }
    throw std::invalid_argument("find_MODS_for_minimum_degree_3 arguments are corrupted");
}

std::vector<NetworKit::node> FominKratschWoegingerDominatingSet::move_to_solution(
        NetworKit::Graph &G, NetworKit::node vertex) {
    std::vector<NetworKit::node> moved_vertices;
    for (auto u : G.neighborRange(vertex)) {
        if (bound.contains(u)) {
            bound.erase(u), free.insert(u);
            moved_vertices.push_back(u);
        }
    }
    return moved_vertices;
}

void FominKratschWoegingerDominatingSet::remove_from_solution(
        NetworKit::node vertex, const std::vector<NetworKit::node> &moved_vertices) {
    for (auto u : std::views::reverse(moved_vertices)) {
        free.erase(u), bound.insert(u);
    }
}

bool FominKratschWoegingerDominatingSet::forget_vertex(
        NetworKit::Graph &G, NetworKit::node vertex, bool is_required) {
    bool is_free = free.contains(vertex);
    if (is_free) {
        free.erase(vertex);
    } else {
        bound.erase(vertex);
    }
    if (G.degree(vertex) == 1) {
        degree[0].erase(vertex);
    } else if (G.degree(vertex) == 2) {
        degree[1].erase(vertex);
    }
    for (auto u : G.neighborRange(vertex)) {
        G.removeEdge(u, vertex);
        if (G.degree(u) == 0) {
            degree[0].erase(u);
        } else if (G.degree(u) == 1) {
            degree[0].insert(u), degree[1].erase(u);
        } else if (G.degree(u) == 2) {
            degree[1].insert(u);
        }
    }
    if (is_required) {
        required.insert(vertex);
    }
    return is_free;
}

void FominKratschWoegingerDominatingSet::retrieve_vertex(
        NetworKit::Graph &G, NetworKit::node vertex, bool is_free, bool is_required) {
    for (auto u : G.neighborRange(vertex)) {
        G.addEdge(u, vertex);
        if (G.degree(u) == 1) {
            degree[0].insert(u);
        } else if (G.degree(u) == 2) {
            degree[0].erase(u), degree[1].insert(u);
        } else if (G.degree(u) == 3) {
            degree[1].erase(u);
        }
    }
    if (G.degree(vertex) == 1) {
        degree[0].insert(vertex);
    } else if (G.degree(vertex) == 2) {
        degree[1].insert(vertex);
    }
    if (is_required) {
        required.erase(vertex);
    }
    if (is_free) {
        free.insert(vertex);
    } else {
        bound.insert(vertex);
    }
}

void SchiermeyerDominatingSet::run() {
    hasRun = true;
    bound.insert(graph->nodeRange().begin(), graph->nodeRange().end());
    NetworKit::Graph core_graph = get_core_graph(*graph, free, bound, required);
    if (bound.empty()) {
        dominating_set.insert(required.begin(), required.end());
        return;
    }
    auto possibilities = merge(free, bound);
    if (find_small_MODS(core_graph, possibilities)) {
        return;
    }
    find_big_MODS(core_graph, possibilities);
}

NetworKit::Graph SchiermeyerDominatingSet::get_core_graph(
        const NetworKit::Graph &G, std::set<NetworKit::node> &free,
        std::set<NetworKit::node> &bound, std::set<NetworKit::node> &required) {
    NetworKit::Graph G_prim(G);
    bool process = true;
    while (process) {
        process = false;
        // Rule 1
        for (auto u : bound) {
            if (G_prim.degree(u) == 0) {
                required.insert(u), process = true;
            }
        }
        std::erase_if(bound, [&G_prim](auto u) { return G_prim.degree(u) == 0; });
        // Rule 2
        std::set<NetworKit::node> rule2change;
        for (auto u : bound) {
            auto it = std::find_if(
                G_prim.neighborRange(u).begin(), G_prim.neighborRange(u).end(),
                [&required](auto v) { return required.contains(v); });
            if (it != G_prim.neighborRange(u).end()) {
                rule2change.insert(u), process = true;
            }
        }
        free.insert(rule2change.begin(), rule2change.end());
        std::erase_if(bound, [&rule2change](auto u) { return rule2change.contains(u); });
        // Rule 3
        std::vector<NetworKit::Edge> edges(G_prim.edgeRange().begin(), G_prim.edgeRange().end());
        for (const auto &[u, v] : edges) {
            if (required.contains(u) || required.contains(v)) {
                G_prim.removeEdge(u, v), process = true;
            } else if (free.contains(u) && free.contains(v)) {
                G_prim.removeEdge(u, v), process = true;
            }
        }
        // Rule 4
        std::set<NetworKit::node> rule4change;
        for (const auto &u : free) {
            auto it = std::find_if(
                G_prim.neighborRange(u).begin(), G_prim.neighborRange(u).end(),
                [&bound](auto v) { return bound.contains(v); });
            if (it != G_prim.neighborRange(u).end()) {
                it = std::find_if(
                    std::next(it), G_prim.neighborRange(u).end(),
                    [&bound](auto v) { return bound.contains(v); });
            }
            if (it == G_prim.neighborRange(u).end()) {
                rule4change.insert(u);
            }
        }
        for (const auto &u : rule4change) {
            G_prim.removeNode(u), G_prim.restoreNode(u), process = true;
        }
        std::erase_if(free, [&rule4change](auto u) { return rule4change.contains(u); });
        // Rule 5
        for (const auto &u : bound) {
            if (G_prim.degree(u) == 1 && !required.contains(u)) {
                auto unique = G_prim.getIthNeighbor(u, 0);
                free.erase(unique), required.insert(unique), process = true;
            }
        }
        std::erase_if(bound, [&required](auto u) { return required.contains(u); });
    }
    return G_prim;
}

bool SchiermeyerDominatingSet::find_small_MODS(
        const NetworKit::Graph &G, const std::vector<NetworKit::node> &V) {
    for (NetworKit::count i = 1; 3 * i <= V.size(); i++) {
        if (find_small_MODS_recursive(G, V, 0, i, dominating_set)) {
            dominating_set.insert(required.begin(), required.end());
            return true;
        }
    }
    return false;
}

void SchiermeyerDominatingSet::find_big_MODS(
        const NetworKit::Graph &G, const std::vector<NetworKit::node> &V) {
    auto solution = find_big_MODS_recursive(G, V, 0);
    dominating_set.insert(solution.begin(), solution.end());
    dominating_set.insert(required.begin(), required.end());
}

std::vector<NetworKit::node> SchiermeyerDominatingSet::find_big_MODS_recursive(
        const NetworKit::Graph &G, const std::vector<NetworKit::node> &V, NetworKit::index index) {
    if (index == V.size()) {
        if (neighborhood.size() < 3 * S.size() || neighborhood.size() >= 3 * (S.size() + 1)) {
            return V;
        }
        for (const auto &u : V) {
            if (S.contains(u)) {
                continue;
            }
            if (neighborhood.size() + get_new_neighborhood(G, u).size() >= 3 * (S.size() + 1)) {
                return V;
            }
        }
        return get_matching_MODS(G, free, bound, S);
    }
    auto solution = find_big_MODS_recursive(G, V, index + 1);
    if (V.size() < 3 * (S.size() + 1)) {
        return solution;
    }
    NetworKit::node next = V.at(index);
    auto added = get_new_neighborhood(G, next);
    S.insert(next), neighborhood.insert(added.begin(), added.end());
    auto other_solution = find_big_MODS_recursive(G, V, index + 1);
    S.erase(next);
    for (auto u : added) {
        neighborhood.erase(u);
    }
    return solution.size() < other_solution.size() ? solution : other_solution;
}

std::vector<NetworKit::node> SchiermeyerDominatingSet::get_new_neighborhood(
        const NetworKit::Graph &G, NetworKit::node vertex) {
    std::vector<NetworKit::node> out;
    if (!neighborhood.contains(vertex)) {
        out.push_back(vertex);
    }
    std::copy_if(
        G.neighborRange(vertex).begin(), G.neighborRange(vertex).end(), std::back_inserter(out),
        [this](auto v){ return !neighborhood.contains(v); });
    return out;
}

std::vector<NetworKit::node> SchiermeyerDominatingSet::get_matching_MODS(
        const NetworKit::Graph &G, std::set<NetworKit::node> free,
        std::set<NetworKit::node> bound, std::set<NetworKit::node> required) {
    auto core_graph = get_core_graph(G, free, bound, required);

    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> boost_graph_t;
    boost_graph_t H(G.numberOfNodes());
    std::map<std::tuple<NetworKit::node, NetworKit::node>, NetworKit::node> owners;
    for (const auto &u : bound) {
        NetworKit::node v = core_graph.getIthNeighbor(u, 0);
        if (v != NetworKit::none && u < v) {
            boost::add_edge(u, v, H);
            owners.emplace(std::make_tuple(u, v), u);
        }
    }
    for (const auto &u : free) {
        NetworKit::node v1 = core_graph.getIthNeighbor(u, 0), v2 = core_graph.getIthNeighbor(u, 1);
        assert(v1 != NetworKit::none && v2 != NetworKit::none);
        if (v1 > v2) {
            std::swap(v1, v2);
        }
        if (!owners.contains(std::make_tuple(v1, v2))) {
            boost::add_edge(v1, v2, H);
            owners.emplace(std::make_tuple(v1, v2), u);
        }
    }
    std::vector<boost::graph_traits<boost_graph_t>::vertex_descriptor> mate(G.numberOfNodes());
    assert(boost::checked_edmonds_maximum_cardinality_matching(H, &mate[0]));

    std::vector<NetworKit::node> optional_dominating_set;
    for (const auto &u : bound) {
        if (mate[u] == NetworKit::none) {
            optional_dominating_set.push_back(u);
        } else if (u < mate[u]) {
            optional_dominating_set.push_back(owners.at(std::make_tuple(u, mate[u])));
        }
    }
    std::copy(required.begin(), required.end(), std::back_inserter(optional_dominating_set));
    return optional_dominating_set;
}

}  /* namespace Koala */
