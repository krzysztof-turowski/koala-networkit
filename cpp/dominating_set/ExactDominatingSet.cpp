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

bool sized_choice_recursive(
        const auto &verify, const std::vector<NetworKit::node> &vertices,
        NetworKit::node index, int remaining, std::vector<bool> &solution) {
    if (index == vertices.size()) {
        return verify(solution);
    }
    if (remaining > 0) {
        solution[vertices.at(index)] = true;
        if (sized_choice_recursive(verify, vertices, index + 1, remaining - 1, solution)) {
            return true;
        }
        solution[vertices.at(index)] = false;
    }
    return index + remaining < vertices.size() && sized_choice_recursive(
        verify, vertices, index + 1, remaining, solution);
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
    dominating_set = recurse();
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
    std::vector<NetworKit::node> possibilities = merge(free, bound);
    for (size_t i = 1; 8 * i <= 3 * possibilities.size(); i++) {
        std::vector<bool> solution(G.numberOfNodes());
        bool found = sized_choice_recursive(
            [&G, this](const std::vector<bool> &S) {
                return is_optional_dominating_set(G, S, bound);
            }, possibilities, 0, i, solution);
        if (found) {
            return solution;
        }
    }
    throw std::invalid_argument("find_MODS_for_minimum_degree_3 arguments are corrupted");
}

std::vector<bool> FominKratschWoegingerDominatingSet::recurse() {
    std::vector<bool> solution;
    if (!degree_one.empty()) {
        NetworKit::node u = *degree_one.begin(), unique = *neighborhood.at(u).begin();
        bool is_free = forget_vertex(u);
        if (is_free) {
            recurse().swap(solution);
        } else {
            auto [moved1, is_free1] = add_to_solution(unique);
            solution = recurse();
            remove_from_solution(unique, moved1, is_free1);
        }
        retrieve_vertex(u, is_free);
        return solution;
    }
    if (!degree_two.empty()) {
        NetworKit::node v = *degree_two.begin();
        std::set<NetworKit::node>::iterator it = neighborhood.at(v).begin();
        NetworKit::node u1 = *it, u2 = *(++it);

        bool is_free = free.contains(v), is1 = forget_vertex(v);
        auto [moved1, is_free1] = add_to_solution(u1);
        recurse().swap(solution);
        remove_from_solution(u1, moved1, is_free1);
        retrieve_vertex(v, is1);
        int solution_size = std::count(solution.begin(), solution.end(), true);

        bool is2 = forget_vertex(u1), is3 = forget_vertex(u2);
        auto [moved2, is_free2] = add_to_solution(v);
        auto solution2 = recurse();
        remove_from_solution(v, moved2, is_free2);
        retrieve_vertex(u2, is3), retrieve_vertex(u1, is2);
        int solution2_size = std::count(solution2.begin(), solution2.end(), true);
        if (solution2_size < solution_size) {
            solution.swap(solution2), solution_size = solution2_size;
        }

        std::vector<bool> solution3;
        bool is4 = forget_vertex(v);
        if (is_free) {
            solution3 = recurse();
        } else {
            auto [moved3, is_free3] = add_to_solution(u2);
            solution3 = recurse();
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
    graph->forNodes([this](NetworKit::node u) {
        bound.insert(u);
    });
    NetworKit::Graph core_graph = core(*graph, free, bound, required);
    if (bound.empty()) {
        graph->forNodes([this](NetworKit::node u) {
            dominating_set.at(u) = required.contains(u);
        });
    }
    possibilities = merge(free, bound);
    if (find_small_MODS(core_graph)) {
        return;
    }
    find_big_MODS(core_graph);
}

bool SchiermeyerDominatingSet::find_small_MODS(const NetworKit::Graph &G) {
    for (int i = 1; 3 * i <= possibilities.size(); i++) {
        bool found = sized_choice_recursive(
            [&G, this](const std::vector<bool> &S) {
                return is_optional_dominating_set(G, S, bound);
            }, possibilities, 0, i, dominating_set);
        if (found) {
            G.forNodes([this](NetworKit::node u) {
                dominating_set.at(u) = dominating_set.at(u) | required.contains(u);
            });
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
    recurse(G, 0);
    for (auto u : required) {
        dominating_set.at(u) = true;
    }
}

void SchiermeyerDominatingSet::recurse(const NetworKit::Graph &G, NetworKit::node index) {
    std::set<NetworKit::node> added;
    if (index == possibilities.size()) {
        if (neighborhood.size() < 3 * choices.size() || neighborhood.size() >= 3 * (choices.size() + 1)) {
            return;
        }
        for (auto u : possibilities) {
            if (choices.contains(u)) {
                continue;
            }
            if (!neighborhood.contains(u)) {
                added.insert(u);
            }
            G.forNeighborsOf(u, [&added, this](NetworKit::node v) {
                if (!neighborhood.contains(v)) {
                    added.insert(v);
                }});
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
    recurse(G, index + 1);

    if (3 * (choices.size() + 1) > possibilities.size()) {
        return;
    }
    NetworKit::node next = possibilities.at(index);
    if (!neighborhood.contains(next)) {
        added.insert(next);
    }
    G.forNeighborsOf(
        next,
        [&added, this] (NetworKit::node v) {
            if (!neighborhood.contains(v)) {
                added.insert(v);
            }
        });

    choices.insert(next);
    for (auto e : added) {
        neighborhood.insert(e);
    }
    recurse(G, index + 1);
    choices.erase(next);
    for (auto e : added) {
        neighborhood.erase(e);
    }
}

std::vector<bool> SchiermeyerDominatingSet::get_matching_MODS(
        const NetworKit::Graph &G, const std::set<NetworKit::node> &required) {
    auto S_prim(required);
    auto T_prim(free);
    auto U_prim(bound);
    auto core_graph = core(G, T_prim, U_prim, S_prim);

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

NetworKit::Graph SchiermeyerDominatingSet::core(
        const NetworKit::Graph &G,
        std::set<NetworKit::node> &free,
        std::set<NetworKit::node> &bound,
        std::set<NetworKit::node> &required) {
    std::vector<std::set<NetworKit::node>> intermediate(G.numberOfNodes());
    G.forEdges([&intermediate](NetworKit::node u, NetworKit::node v) {
        intermediate[u].insert(v);
        intermediate[v].insert(u);
    });
    bool process = true;
    while (process) {
        process = false;
        for (auto e : bound) {
            if (intermediate[e].empty()) {
                process = true;
                required.insert(e);
            }
        }
        std::erase_if(
            bound,
            [&intermediate](NetworKit::node u) {
                return intermediate[u].empty();
            });
        std::set<NetworKit::node> rule2change;
        for (auto e : bound) {
            for (auto nei : intermediate[e]) {
                if (required.contains(nei)) {
                    process = true;
                    rule2change.insert(e);
                    break;
                }
            }
        }
        for (auto e : rule2change) {
            free.insert(e);
        }
        std::erase_if(bound, [&rule2change](NetworKit::node u) {return rule2change.contains(u);});
        for (int i = 0; i < intermediate.size(); i++) {
            if (free.contains(i)) {
                auto erased = std::erase_if(
                    intermediate[i],
                    [&free, &required](NetworKit::node e) {
                        return free.contains(e) || required.contains(e);
                    });
                if (erased > 0) {
                    process = true;
                }
            } else if (required.contains(i)) {
                if (!intermediate[i].empty()) {
                    process = true;
                    for (auto e : intermediate[i]) {
                        if (bound.contains(e)) {
                            intermediate.at(e).erase(i);
                        }
                    }
                    intermediate[i].clear();
                }
            }
        }
        std::set<NetworKit::node> rule4change;
        for (auto e : free) {
            int boundCount = 0;
            for (auto nei : intermediate[e]) {
                if (bound.contains(nei)) {
                    boundCount++;
                }
                if (boundCount > 1) {
                    break;
                }
            }
            if (boundCount <= 1) {
                process = true;
                rule4change.insert(e);
            }
        }
        for (auto e : rule4change) {
            for (auto nei : intermediate[e]) {
                intermediate[nei].erase(e);
            }
            intermediate[e].clear();
        }
        std::erase_if(free, [&rule4change](NetworKit::node u) { return rule4change.contains(u); });
        for (auto e : bound) {
            if (intermediate[e].size() == 1 && !required.contains(e)) {
                process = true;
                auto unique = *intermediate[e].begin();
                free.erase(unique);
                required.insert(unique);
            }
        }
        std::erase_if(bound, [&required](NetworKit::node u) { return required.contains(u); });
    }

    NetworKit::Graph graph(intermediate.size());
    for (NetworKit::node i = 0; i < intermediate.size(); i++) {
        for (auto e : intermediate.at(i)) {
            assert(intermediate.at(e).contains(i));
            if (i < e) {
                graph.addEdge(i, e);
            }
        }
    }
    return graph;
}

}  /* namespace Koala */
