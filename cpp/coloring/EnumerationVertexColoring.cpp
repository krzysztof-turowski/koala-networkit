/*
 * EnumerationVertexColoring.cpp
 *
 * Created on: 06.02.2023
 *   Author: Zofia Glapa (zofia.glapa@student.uj.edu.pl)
 */

#include <coloring/ExactVertexColoring.hpp>
#include <iostream>
#include <limits>
#include <networkit/auxiliary/BucketPQ.hpp>

namespace Koala {

EnumerationVertexColoring::EnumerationVertexColoring(const NetworKit::Graph& graph)
    : graph(std::make_optional(graph)) {
}

const std::map<NetworKit::node, int> EnumerationVertexColoring::getColoring() const {
    assureFinished();
    std::map<NetworKit::node, int> final_coloring;
    if (final_coloring.empty()) {
        for (int i = 0; i < graph->numberOfNodes(); ++i) {
            final_coloring[ordering[i]] = best_solution[i];
        }
    }
    return final_coloring;
}

void EnumerationVertexColoring::forwards() {
    for (int i = r; i < graph->numberOfNodes(); ++i) {
        if (r == 0 || i >= r + 1) {
            determine_feasible_colors(i);
        }
        if (feasible_colors[i].empty()) {
            r = i;
            return;
        }
        current_solution[i] = *(feasible_colors[i]).begin();
    }
    best_solution = current_solution;
    NetworKit::node maximal_color_index = NetworKit::none;
    int maximal_color = 0;
    for (int i = 0; i < graph->numberOfNodes(); ++i) {
        if (current_solution[i] > maximal_color) {            
            maximal_color_index = i, maximal_color = current_solution[i];
        }
    }
    r = maximal_color_index, current_bound = maximal_color;
}

void EnumerationVertexColoring::backwards() {
    determine_current_predecessors(r);
    while (!current_predecessors.empty()) {
        int i = *current_predecessors.begin();
        current_predecessors.erase(i);
        feasible_colors[i].erase(current_solution[i]);
        if (!feasible_colors[i].empty() && *(feasible_colors[i]).begin() < current_bound) {
            r = i;
            return;
        }
    }
    r = 0;
}

void EnumerationVertexColoring::determine_feasible_colors(int i) {
    std::set<int> feasible_colors_for_node;
    int current_maximal_color = 0;
    for (int j = 0; j < i; ++j) {
        if (current_maximal_color < current_solution[j]) {
            current_maximal_color = current_solution[j];
        }
    }
    for (int j = 1; j <= current_maximal_color + 1; ++j) {
        if (j >= current_bound)
            break;
        feasible_colors_for_node.insert(j);
    }
    for (int j = 0; j < i; ++j) {
        if (graph->hasEdge(ordering[i], ordering[j])) {
            feasible_colors_for_node.erase(current_solution[j]);
        }
    }
    feasible_colors[i] = feasible_colors_for_node;
}

std::vector<NetworKit::node> BrownEnumerationVertexColoring::greedy_largest_first_ordering() {
    std::unordered_set<NetworKit::node> already_ordered;
    std::vector<NetworKit::node> ordering;
    std::map<NetworKit::node, int> number_of_neighbours_in_ordering;
    auto compare = [](const std::tuple<int, int, NetworKit::node>& a,
                   const std::tuple<int, int, NetworKit::node>& b) {
        if (std::get<0>(a) == std::get<0>(b)) {
            if (std::get<1>(a) == std::get<1>(b)) {
                return std::get<2>(a) < std::get<2>(b);
            }
            return std::get<1>(a) > std::get<1>(b);
        }
        return std::get<0>(a) > std::get<0>(b);
    };
    std::set<std::tuple<int, int, NetworKit::node>, decltype(compare)> number_of_neighbours_in_ordering_queue;
    NetworKit::node current_node;

    NetworKit::node max_degree_node = [](const NetworKit::Graph& graph) {
        NetworKit::node max_degree_node = 0;
        unsigned int max_degree = 0;
        graph.forNodes([&](NetworKit::node u) {
            if (graph.degree(u) > max_degree) {
                max_degree = graph.degree(u);
                max_degree_node = u;
            }
        });
        return max_degree_node;
    }(*graph);

    ordering.push_back(max_degree_node);
    already_ordered.insert(max_degree_node);

    graph->forNodes([&](NetworKit::node u) {
        if (u == max_degree_node) {
            return;
        }
        number_of_neighbours_in_ordering[u] = graph->hasEdge(u, max_degree_node);
        number_of_neighbours_in_ordering_queue.insert(
        std::make_tuple(number_of_neighbours_in_ordering[u], graph->degree(u), u));
    });
    while (already_ordered.size() < graph->numberOfNodes()) {
        current_node = get<2>(*number_of_neighbours_in_ordering_queue.begin());
        number_of_neighbours_in_ordering_queue.erase(
        number_of_neighbours_in_ordering_queue.begin());

        ordering.push_back(current_node);
        already_ordered.insert(current_node);

        graph->forNeighborsOf(current_node, [&](NetworKit::node v) {
            if (already_ordered.find(v) == already_ordered.end()) {
                number_of_neighbours_in_ordering[v]++;
                number_of_neighbours_in_ordering_queue.erase(
                std::make_tuple(number_of_neighbours_in_ordering[v] - 1, graph->degree(v), v));
                number_of_neighbours_in_ordering_queue.insert(
                std::make_tuple(number_of_neighbours_in_ordering[v], graph->degree(v), v));
            }
        });
    }
    return ordering;
}

void BrownEnumerationVertexColoring::determine_current_predecessors(int r) {
    current_predecessors = std::set<int, std::greater<int>>();
    for (int i = 0; i < r; ++i) {
        current_predecessors.insert(i);
    }
}

void BrownEnumerationVertexColoring::run() {
    ordering = greedy_largest_first_ordering();
    lower_bound = 1, upper_bound = graph->numberOfNodes();
    r = 0;
    current_bound = upper_bound + 1;
    feasible_colors.resize(graph->numberOfNodes());
    feasible_colors[0].insert(1);
    current_solution.resize(graph->numberOfNodes());

    while (true) {
        forwards();
        if (current_bound == lower_bound) {
            break;
        }
        backwards();
        if (r == 0) {
            break;
        }
    }
    hasRun = true;
}

void ChristofidesEnumerationVertexColoring::calculate_transitive_closure() {
    // use APSP
    transitive_closure.resize(graph->numberOfNodes());
    for (int u = 0; u < graph->numberOfNodes(); ++u) {
        transitive_closure[u].resize(graph->numberOfNodes());
    }
    for (int u = 0; u < graph->numberOfNodes(); ++u) {
        for (int v = 0; v < graph->numberOfNodes(); ++v) {
            if (graph->hasEdge(ordering[u], ordering[v]) && u < v) {
                transitive_closure[u][v] = true;
            }
        }
    }
    for (int u = 0; u < graph->numberOfNodes(); ++u) {
        for (int v = 0; v < graph->numberOfNodes(); ++v) {
            if (transitive_closure[u][v]) {
                for (int w = 0; w < graph->numberOfNodes(); ++w) {
                    if (transitive_closure[v][w]) {
                        transitive_closure[u][w] = true;
                    }
                }
            }
        }
    }
}

void ChristofidesEnumerationVertexColoring::determine_current_predecessors(int r) {
    for (int u = 0; u < graph->numberOfNodes(); ++u) {
        if (transitive_closure[u][r]) {
            current_predecessors.insert(u);
        }
    }
}

const std::vector<std::vector<bool>>&
ChristofidesEnumerationVertexColoring::getTransitiveClosure() const {
    assureFinished();
    return transitive_closure;
}

void ChristofidesEnumerationVertexColoring::run() {
    ordering = greedy_largest_first_ordering();
    lower_bound = 1, upper_bound = graph->numberOfNodes();
    r = 0;
    current_bound = upper_bound + 1;

    calculate_transitive_closure();

    feasible_colors.resize(graph->numberOfNodes());
    feasible_colors[0].insert(1);
    current_solution.resize(graph->numberOfNodes());

    while (true) {
        forwards();
        if (current_bound == lower_bound) {
            break;
        }
        backwards();
        if (r == 0) {
            break;
        }
    }

    hasRun = true;
}

std::vector<NetworKit::node> BrelazEnumerationVertexColoring::interchange_component(
std::vector<NetworKit::node>& scurrent_boundgraph,
std::map<NetworKit::node, int>& solution,
NetworKit::node new_node,
int alpha) {
    std::unordered_set<NetworKit::node> visited;
    std::vector<NetworKit::node> verticesToRecolor;
    for (NetworKit::node u : scurrent_boundgraph) {
        if (visited.find(u) == visited.end()) {
            std::vector<NetworKit::node> component;
            std::queue<NetworKit::node> queue;
            std::unordered_set<int> neighborsColors;
            queue.push(u);
            visited.insert(u);
            while (!queue.empty()) {
                NetworKit::node v = queue.front();
                queue.pop();
                if (graph->hasEdge(new_node, v)) {
                    if (neighborsColors.size() == 0) {
                        neighborsColors.insert(solution[v]);
                    } else if (neighborsColors.find(solution[v]) == neighborsColors.end()) {
                        return std::vector<NetworKit::node>();
                    }
                }
                component.push_back(v);
                for (NetworKit::node w : scurrent_boundgraph) {
                    if (graph->hasEdge(v, w) && visited.find(w) == visited.end()) {
                        queue.push(w);
                        visited.insert(w);
                    }
                }
            }
            if (neighborsColors.find(alpha) != neighborsColors.end()) {
                verticesToRecolor.insert(
                verticesToRecolor.end(), component.begin(), component.end());
            }
        }
    }
    return verticesToRecolor;
}

bool BrelazEnumerationVertexColoring::is_interchangeable(
        std::vector<int>& K, NetworKit::node new_node, std::map<NetworKit::node, int>& solution) {
    for (int alpha : K) {
        for (int beta : K) {
            if (alpha == beta) {
                break;
            }
            std::vector<NetworKit::node> scurrent_boundgraph;
            for (const auto& [node, color] : solution) {
                if (color == alpha || color == beta) {
                    scurrent_boundgraph.push_back(node);
                }
            }
            auto interchangeable_component =
                interchange_component(scurrent_boundgraph, solution, new_node, alpha);
            if (interchangeable_component.size() > 0) {
                for (NetworKit::node v : interchangeable_component) {
                    solution[v] = (alpha == solution[v]) ? beta : alpha;
                }
                solution[new_node] = alpha;
                return true;
            }
        }
    }
    return false;
}

std::vector<NetworKit::node> BrelazEnumerationVertexColoring::saturation_largest_first_with_interchange() {
    std::map<NetworKit::node, int> solution;
    std::vector<NetworKit::node> ordering;
    NetworKit::node max_degree_node = 0;
    uint max_degree = 0;

    graph->forNodes([&](NetworKit::node u) {
        if (graph->degree(u) > max_degree) {
            max_degree = graph->degree(u);
            max_degree_node = u;
        }
    });

    solution.insert(std::make_pair(max_degree_node, 1));
    ordering.push_back(max_degree_node);

    auto satur_comp = [&](const auto &a, const auto &b) {
        if (std::get<0>(a) == std::get<0>(b)) {
            if (std::get<1>(a) == std::get<1>(b)) {
                return std::get<2>(a) < std::get<2>(b);
            }
            return std::get<1>(a) > std::get<1>(b);
        }
        return std::get<0>(a) > std::get<0>(b);
    };

    auto saturation =
        std::set<std::tuple<int, int, NetworKit::node>, decltype(satur_comp)>(satur_comp);

    std::unordered_map<NetworKit::node, std::unordered_set<int>> neighbours_colors;
    graph->forNodes([&](NetworKit::node u) {
        if (u != max_degree_node) {
            if (graph->hasEdge(max_degree_node, u)) {
                saturation.insert(std::make_tuple(1, graph->degree(u) - 1, u));
                neighbours_colors[u].insert(1);
            } else {
                saturation.insert(std::make_tuple(0, graph->degree(u), u));
            }
        }
    });

    int max_color = 1, max_clique_size = 0;
    while (solution.size() < graph->numberOfNodes()) {
        auto max_saturation = saturation.begin();
        NetworKit::node u = std::get<2>(*max_saturation);
        saturation.erase(max_saturation);
        ordering.push_back(u);

        int first_valid_color = 1;
        std::set<int> forbidden_colors;
        std::vector<std::pair<int, std::vector<NetworKit::node>>> neighbours_color_count;
        neighbours_color_count.resize(max_color + 1);
        graph->forNeighborsOf(u, [&](NetworKit::node v) {
            if (solution.find(v) != solution.end()) {
                forbidden_colors.insert(solution[v]);
                if (neighbours_color_count[solution[v]].first == 0) {
                    neighbours_color_count[solution[v]].first = 1;
                    neighbours_color_count[solution[v]].second.push_back(v);
                }
            }
        });
        while (forbidden_colors.find(first_valid_color) != forbidden_colors.end()) {
            ++first_valid_color;
        }

        if (first_valid_color <= max_color) {
            solution.insert(std::make_pair(u, first_valid_color));
            if (max_clique_size == 0) {
                max_clique_size = max_color;
            }
        } else {
            std::vector<int> K;
            for (int i = 1; i <= max_color; ++i) {
                if (neighbours_color_count[i].first == 1) {
                    K.push_back(i);
                }
            }
            auto interchange = is_interchangeable(K, u, solution);
            if (!interchange) {
                solution.insert(std::make_pair(u, ++max_color));
            } else {
                if (max_clique_size == 0) {
                    max_clique_size = max_color;
                }
            }
        }
        std::vector<std::tuple<int, int, NetworKit::node>> temp_satur;
        for (const auto& [saturation_degree, degree, node] : saturation) {
            if (graph->hasEdge(u, node)) {
                neighbours_colors[node].insert(solution[u]);
                temp_satur.push_back(std::make_tuple(saturation_degree, degree, node));
            }
        }
        for (const auto& [saturation_degree, degree, node] : temp_satur) {
            saturation.erase(std::make_tuple(saturation_degree, degree, node));
            saturation.insert(std::make_tuple(neighbours_colors[node].size(), degree - 1, node));
        }
    }
    lower_bound = max_clique_size, upper_bound = max_color;

    for (int i = 0; i < graph->numberOfNodes(); ++i) {
        current_solution[i] = solution[ordering[i]];
    }
    best_solution = current_solution;
    return ordering;
}

std::vector<int> BrelazEnumerationVertexColoring::get_representatives_of_adjacent_predecessors(
        int i) {
    std::vector<int> representatives(current_bound, INT_MAX);
    for (int j = 0; j < i; ++j) {
        if (graph->hasEdge(ordering[j], ordering[i]) && current_solution[j] < current_bound) {
            if (representatives[current_solution[j]] > j) {
                representatives[current_solution[j]] = j;
            }
        }
    }
    return representatives;
}

void BrelazEnumerationVertexColoring::determine_current_predecessors(int r) {
    auto representatives = get_representatives_of_adjacent_predecessors(r);
    for (int representative : representatives) {
        if (representative < INT_MAX) {
            current_predecessors.insert(representative);
        }
    }
}

void BrelazEnumerationVertexColoring::backwards() {
    determine_current_predecessors(r);
    while (!current_predecessors.empty()) {
        int i = *current_predecessors.begin();
        current_predecessors.erase(i);
        determine_current_predecessors(i);
        feasible_colors[i].erase(current_solution[i]);
        if (!feasible_colors[i].empty()) {
            r = i;
            return;
        }
    }
    r = 0;
}

void BrelazEnumerationVertexColoring::run() {
    current_solution.resize(graph->numberOfNodes());
    best_solution.resize(graph->numberOfNodes());

    ordering = saturation_largest_first_with_interchange();

    if (lower_bound != upper_bound) {
        r = 0;
        current_bound = upper_bound;

        feasible_colors.resize(graph->numberOfNodes());
        feasible_colors[0].insert(1);

        while (true) {
            forwards();
            if (current_bound == lower_bound) {
                break;
            }
            backwards();
            if (r == 0) {
                break;
            }
        }
    }
    hasRun = true;
}

void KormanEnumerationVertexColoring::forwards() {
    std::vector<bool> is_colored(graph->numberOfNodes(), false);
    std::vector<std::unordered_set<int>> neighbour_colors(graph->numberOfNodes());
    Aux::BucketPQ queue(graph->numberOfNodes(), -graph->numberOfNodes() + 1, 0);

    for (int i : new_ordering) {
        is_colored[i] = true;
        graph->forNeighborsOf(ordering[i], [&](NetworKit::node v) {
            auto j = position[v];
            neighbour_colors[j].insert(current_solution[i]);
        });
    }

    graph->forNodes([&](NetworKit::node u) {
        auto j = position[u];
        if (!is_colored[j]) {
            queue.insert(-neighbour_colors[j].size(), j);
        }
    });

    while (!queue.empty()) {
        auto top = queue.extractMin();
        auto node = top.second;

        new_ordering.push_back(node);

        determine_feasible_colors(new_ordering.size() - 1, neighbour_colors[node]);

        if (feasible_colors[node].empty()) {
            r = new_ordering.size() - 1;
            return;
        }

        is_colored[node] = true;
        current_solution[node] = *feasible_colors[node].begin();

        graph->forNeighborsOf(ordering[node], [&](NetworKit::node v) {
            auto j = position[v];
            if (!is_colored[j]) {
                neighbour_colors[j].insert(current_solution[node]);
                queue.changeKey(-neighbour_colors[j].size(), j);
            }
        });
    }

    best_solution = current_solution;
    int maximal_color = 0;
    int maximal_color_index = -1;
    for (int i = 0; i < graph->numberOfNodes(); ++i) {
        if (best_solution[new_ordering[i]] > maximal_color) {
            maximal_color = best_solution[new_ordering[i]];
            maximal_color_index = i;
        }
    }
    current_bound = maximal_color;
    r = maximal_color_index;
}

void KormanEnumerationVertexColoring::backwards() {
    for (auto i = r - 1; i >= 0; i--) {
        feasible_colors[new_ordering[i]].erase(current_solution[new_ordering[i]]);
        if (!feasible_colors[new_ordering[i]].empty()) {
            if (*feasible_colors[new_ordering[i]].begin() < current_bound) {
                current_solution[new_ordering[i]] = *feasible_colors[new_ordering[i]].begin();
                while (new_ordering.size() > i + 1)
                    new_ordering.pop_back();
                r = i;
                return;
            }
        }
    }
    r = 0;
}

void KormanEnumerationVertexColoring::determine_feasible_colors(int i,
std::unordered_set<int> blocked_colors) {
    std::set<int> feasible_colors_for_node;
    int current_maximal_color = 0;

    for (int j = 0; j < i; ++j) {
        if (current_maximal_color < current_solution[new_ordering[j]]) {
            current_maximal_color = current_solution[new_ordering[j]];
        }
    }

    for (int j = 1; j <= current_maximal_color + 1; ++j) {
        if (j >= current_bound)
            break;
        if (blocked_colors.find(j) == blocked_colors.end()) {
            feasible_colors_for_node.insert(j);
        }
    }

    feasible_colors[new_ordering[i]] = feasible_colors_for_node;
}

void KormanEnumerationVertexColoring::run() {
    ordering = greedy_largest_first_ordering();
    lower_bound = 1;
    upper_bound = graph->numberOfNodes();
    current_solution.resize(graph->numberOfNodes());
    best_solution.resize(graph->numberOfNodes());

    for (int i = 0; i < graph->numberOfNodes(); i++) {
        position[ordering[i]] = i;
    }

    r = 0;
    current_bound = upper_bound + 1;

    new_ordering.push_back(0);
    feasible_colors.resize(graph->numberOfNodes());
    feasible_colors[0].insert(1);
    current_solution[0] = 1;

    while (true) {
        forwards();
        if (current_bound == lower_bound) {
            break;
        }
        backwards();
        if (r == 0) {
            break;
        }
    }
    hasRun = true;
}

} /* namespace Koala */
