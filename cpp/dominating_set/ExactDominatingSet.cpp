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
        const NetworKit::Graph &G, const std::vector<bool> &solution,
        std::set<NetworKit::node> bound) {
    for (int u = 0; u < solution.size(); u++) {
        if (solution[u]) {
            bound.erase(u);
            G.forNeighborsOf(u, [&bound](NetworKit::node v) { bound.erase(v); });
        }
    }
    return bound.empty();
}

std::vector<NetworKit::node> merge(
        const std::set<NetworKit::node> &A, const std::set<NetworKit::node> &B) {
    std::vector<NetworKit::node> merged;
    std::set_union(A.begin(), A.end(), B.begin(), B.end(), std::inserter(merged, merged.begin()));
    return merged;
}

bool find_small_MODS_recursive(
        const NetworKit::Graph &G, const std::set<NetworKit::node> &bound,
        const std::vector<NetworKit::node> &vertices, NetworKit::count index,
        NetworKit::count remaining, std::vector<bool> &solution) {
    if (index == vertices.size()) {
        return is_optional_dominating_set(G, solution, bound);
    }
    if (remaining > 0) {
        solution[vertices.at(index)] = true;
        if (find_small_MODS_recursive(G, bound, vertices, index + 1, remaining - 1, solution)) {
            return true;
        }
        solution[vertices.at(index)] = false;
    }
    return index + remaining < vertices.size() && find_small_MODS_recursive(
        G, bound, vertices, index + 1, remaining, solution);
}

void FominKratschWoegingerDominatingSet::run() {
    hasRun = true;
    neighborhood.resize(graph->numberOfNodes());
    graph->forNodes([this](NetworKit::node u) {
        bound.insert(u);
        graph->forNeighborsOf(u, [u, this](NetworKit::node neighbor) {
            neighborhood[u].insert(neighbor);
        });
        if (graph->degree(u) == 1) {
            degree_one.insert(u);
        } else if (graph->degree(u) == 2) {
            degree_two.insert(u);
        }
    });
    dominating_set = find_big_MODS_recursive();
}

std::vector<bool> FominKratschWoegingerDominatingSet::find_MODS_for_minimum_degree_3() {
    NetworKit::Graph G(neighborhood.size());
    for (NetworKit::node u = 0; u < neighborhood.size(); u++) {
        if (!required.contains(u)) {
          for (NetworKit::node v : neighborhood.at(u)) {
              if (!required.contains(v)) {
                  G.addEdge(u, v);
              }
          }
      }
    }
    auto possibilities = merge(free, bound);
    for (size_t i = 1; 8 * i <= 3 * possibilities.size(); i++) {
        std::vector<bool> solution(G.numberOfNodes());
        if (find_small_MODS_recursive(G, bound, possibilities, 0, i, solution)) {
            return solution;
        }
    }
    throw std::invalid_argument("find_MODS_for_minimum_degree_3 arguments are corrupted");
}

std::vector<bool> FominKratschWoegingerDominatingSet::find_big_MODS_recursive() {
    std::vector<bool> solution;
    if (!degree_one.empty()) {
        NetworKit::node u = *degree_one.begin(), unique = *neighborhood.at(u).begin();
        bool is_free = forget_vertex(u);
        if (is_free) {
            find_big_MODS_recursive().swap(solution);
        } else {
            auto [moved1, is_free1] = add_to_solution(unique);
            solution = find_big_MODS_recursive();
            remove_from_solution(unique, moved1, is_free1);
        }
        retrieve_vertex(u, is_free);
        return solution;
    }
    if (!degree_two.empty()) {
        NetworKit::node v = *degree_two.begin();
        auto it = neighborhood.at(v).begin();
        NetworKit::node u1 = *it, u2 = *(++it);

        bool is_free = free.contains(v), is1 = forget_vertex(v);
        auto [moved1, is_free1] = add_to_solution(u1);
        find_big_MODS_recursive().swap(solution);
        remove_from_solution(u1, moved1, is_free1);
        retrieve_vertex(v, is1);
        int solution_size = std::count(solution.begin(), solution.end(), true);

        bool is2 = forget_vertex(u1), is3 = forget_vertex(u2);
        auto [moved2, is_free2] = add_to_solution(v);
        auto solution2 = find_big_MODS_recursive();
        remove_from_solution(v, moved2, is_free2);
        retrieve_vertex(u2, is3), retrieve_vertex(u1, is2);
        int solution2_size = std::count(solution2.begin(), solution2.end(), true);
        if (solution2_size < solution_size) {
            solution.swap(solution2), solution_size = solution2_size;
        }

        std::vector<bool> solution3;
        bool is4 = forget_vertex(v);
        if (is_free) {
            solution3 = find_big_MODS_recursive();
        } else {
            auto [moved3, is_free3] = add_to_solution(u2);
            solution3 = find_big_MODS_recursive();
            remove_from_solution(u2, moved3, is_free3);
        }
        retrieve_vertex(v, is4);
        int solution3_size = std::count(solution3.begin(), solution3.end(), true);
        if (solution3_size < solution_size) {
            solution.swap(solution3), solution_size = solution3_size;
        }
        return solution;
    }
    std::set<NetworKit::node> bound_isolated_vertices;
    for (int i = 0; i < neighborhood.size(); i++) {
        if (neighborhood.at(i).empty() && bound.contains(i)) {
            bound.erase(i), bound_isolated_vertices.insert(i);
        }
    }
    if (bound.empty()) {
        solution.resize(neighborhood.size());
    } else {
        find_MODS_for_minimum_degree_3().swap(solution);
    }
    for (auto u : required) {
        solution.at(u) = true;
    }
    for (auto u : bound_isolated_vertices) {
        solution.at(u) = true;
        bound.insert(u);
    }
    return solution;
}

std::tuple<std::vector<NetworKit::node>, bool> FominKratschWoegingerDominatingSet::add_to_solution(
        NetworKit::node vertex) {
    std::vector<NetworKit::node> moved_vertices;
    for (auto neighbor : neighborhood.at(vertex)) {
        if (bound.contains(neighbor)) {
            bound.erase(neighbor), free.insert(neighbor);
            moved_vertices.push_back(neighbor);
        }
    }
    bool is_free = forget_vertex(vertex);
    required.insert(vertex);
    return {moved_vertices, is_free};
}

void FominKratschWoegingerDominatingSet::remove_from_solution(
        NetworKit::node vertex, std::vector<NetworKit::node> &moved_vertices, bool is_free) {
    required.erase(vertex);
    retrieve_vertex(vertex, is_free);
    for (auto neighbor : std::views::reverse(moved_vertices)) {
        free.erase(neighbor);
        bound.insert(neighbor);
    }
}

bool FominKratschWoegingerDominatingSet::forget_vertex(NetworKit::node vertex) {
    bool is_free = free.contains(vertex);
    if (is_free) {
        free.erase(vertex);
    } else {
        bound.erase(vertex);
    }
    auto degree = neighborhood.at(vertex).size();
    if (degree == 1) {
        degree_one.erase(vertex);
    } else if (degree == 2) {
        degree_two.erase(vertex);
    }
    for (auto neighbor : neighborhood.at(vertex)) {
        neighborhood.at(neighbor).erase(vertex);
        on_degree_decrement(neighbor);
    }
    return is_free;
}

void FominKratschWoegingerDominatingSet::retrieve_vertex(NetworKit::node vertex, bool is_free) {
    for (auto neighbor : neighborhood.at(vertex)) {
        neighborhood.at(neighbor).insert(vertex);
        on_degree_increment(neighbor);
    }
    auto degree = neighborhood.at(vertex).size();
    if (degree == 1) {
        degree_one.insert(vertex);
    } else if (degree == 2) {
        degree_two.insert(vertex);
    }
    if (is_free) {
        free.insert(vertex);
    } else {
        bound.insert(vertex);
    }
}

void FominKratschWoegingerDominatingSet::on_degree_decrement(NetworKit::node vertex) {
    auto degree = neighborhood.at(vertex).size();
    if (degree == 0) {
        degree_one.erase(vertex);
    } else if (degree == 1) {
        degree_one.insert(vertex), degree_two.erase(vertex);
    } else if (degree == 2) {
        degree_two.insert(vertex);
    }
}

void FominKratschWoegingerDominatingSet::on_degree_increment(NetworKit::node vertex) {
    auto degree = neighborhood.at(vertex).size();
    if (degree == 1) {
        degree_one.insert(vertex);
    } else if (degree == 2) {
        degree_one.erase(vertex), degree_two.insert(vertex);
    } else if (degree == 3) {
        degree_two.erase(vertex);
    }
}

void SchiermeyerDominatingSet::run() {
    hasRun = true;
    dominating_set.resize(graph->numberOfNodes());
    bound.insert(graph->nodeRange().begin(), graph->nodeRange().end());
    NetworKit::Graph core_graph = get_core(*graph, free, bound, required);
    if (bound.empty()) {
        for (auto u : required) {
            dominating_set.at(u) = true;
        }
    }
    possibilities = merge(free, bound);
    if (find_small_MODS(core_graph)) {
        return;
    }
    find_big_MODS(core_graph);
}

bool SchiermeyerDominatingSet::find_small_MODS(const NetworKit::Graph &G) {
    for (int i = 1; 3 * i <= possibilities.size(); i++) {
        if (find_small_MODS_recursive(G, bound, possibilities, 0, i, dominating_set)) {
            for (auto u : required) {
                dominating_set.at(u) = true;
            }
            return true;
        }
    }
    return false;
}

void SchiermeyerDominatingSet::find_big_MODS(const NetworKit::Graph &G) {
    dominating_set = std::vector<bool>(G.numberOfNodes());
    for (auto u : possibilities) {
        dominating_set.at(u) = true;
    }
    find_big_MODS_recursive(G, 0);
    for (auto u : required) {
        dominating_set.at(u) = true;
    }
}

void SchiermeyerDominatingSet::find_big_MODS_recursive(
        const NetworKit::Graph &G, NetworKit::count index) {
    if (index == possibilities.size()) {
        if (neighborhood.size() < 3 * choices.size()) {
            return;
        }
        if (neighborhood.size() >= 3 * (choices.size() + 1)) {
            return;
        }
        for (auto u : possibilities) {
            if (choices.contains(u)) {
                continue;
            }
            auto added = get_new_neighborhood(G, u);
            if (neighborhood.size() + added.size() >= 3 * (choices.size() + 1)) {
                return;
            }
        }
        auto optional_dominating_set = get_matching_MODS(G, choices);
        auto optional_dominating_set_size = std::count(optional_dominating_set.begin(), optional_dominating_set.end(), true);
        auto dominating_set_size = std::count(dominating_set.begin(), dominating_set.end(), true);
        if (optional_dominating_set_size < dominating_set_size) {
            dominating_set = optional_dominating_set;
        }
        return;
    }
    find_big_MODS_recursive(G, index + 1);
    if (3 * (choices.size() + 1) > possibilities.size()) {
        return;
    }
    NetworKit::node next = possibilities.at(index);
    auto added = get_new_neighborhood(G, next);
    choices.insert(next), neighborhood.insert(added.begin(), added.end());
    find_big_MODS_recursive(G, index + 1);
    choices.erase(next);
    for (auto u : added) {
        neighborhood.erase(u);
    }
}

std::set<NetworKit::node> SchiermeyerDominatingSet::get_new_neighborhood(
        const NetworKit::Graph &G, NetworKit::node vertex) {
    std::set<NetworKit::node> out;
    if (!neighborhood.contains(vertex)) {
        out.insert(vertex);
    }
    std::copy_if(
        G.neighborRange(vertex).begin(), G.neighborRange(vertex).end(),
        std::inserter(out, out.begin()),
        [this](auto v){ return !neighborhood.contains(v); });
    return out;
}

std::vector<bool> SchiermeyerDominatingSet::get_matching_MODS(
        const NetworKit::Graph &G, const std::set<NetworKit::node> &required) {
    auto S_prim(required);
    auto T_prim(free);
    auto U_prim(bound);
    auto core_graph = get_core(G, T_prim, U_prim, S_prim);

    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> boost_graph_t;
    boost_graph_t H(G.numberOfNodes());
    std::map<std::tuple<NetworKit::node, NetworKit::node>, NetworKit::node> owners;
    for (auto u : U_prim) {
        NetworKit::node v = core_graph.getIthNeighbor(u, 0);
        if (v != NetworKit::none && u < v) {
            boost::add_edge(u, v, H);
            owners.emplace(std::make_tuple(u, v), u);
        }
    }
    for (auto u : T_prim) {
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

    std::vector<bool> optional_dominating_set(G.numberOfNodes());
    for (auto u : U_prim) {
        if (mate[u] == NetworKit::none) {
            optional_dominating_set[u] = true;
        } else if (u < mate[u]) {
            optional_dominating_set[owners.at(std::make_tuple(u, mate[u]))] = true;
        }
    }
    for (auto e : S_prim) {
        optional_dominating_set[e] = true;
    }
    return optional_dominating_set;
}

NetworKit::Graph SchiermeyerDominatingSet::get_core(
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
                auto unique = *(G_prim.neighborRange(u).begin());
                free.erase(unique), required.insert(unique), process = true;
            }
        }
        std::erase_if(bound, [&required](auto u) { return required.contains(u); });
    }
    return G_prim;
}

}  /* namespace Koala */
