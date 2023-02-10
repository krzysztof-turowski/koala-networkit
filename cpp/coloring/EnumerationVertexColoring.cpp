/*
* EnumerationVertexColoring.cpp
*
* Created on: 06.02.2023
*   Author: Zofia Glapa (zofia.glapa@student.uj.edu.pl)
*/

#include <coloring/EnumerationVertexColoring.hpp>
#include <limits>

namespace Koala {

EnumerationVertexColoring::EnumerationVertexColoring(const NetworKit::Graph &graph) :
        graph(std::make_optional(graph)) {}

const std::map<NetworKit::node, int>& EnumerationVertexColoring::getColoring() const {
    assureFinished();
    return best_solution;
}

void EnumerationVertexColoring::forwards() {
    for (int i = r; i < n; ++i) {
        if (r < i)
            determine_feasible_colors(i);
        if (feasible_colors[ordering[i]].empty()) {
            r = i;
            return;
        }
        current_solution[ordering[i]] = *feasible_colors[ordering[i]].begin();
    }
    best_solution = current_solution;
    int maximal_color = 0;
    NetworKit::node maximal_color_node = -1;
    int maximal_color_index = -1;
    for (auto& node_color : best_solution) {
        if (node_color.second > maximal_color) {
            maximal_color = node_color.second;
            maximal_color_node = node_color.first;
        }
    }

    for (int i = 0; i < n; ++i) {
        if (ordering[i] == maximal_color_node) {
            maximal_color_index = i;
            break;
        }
    }
    ub = maximal_color;
    r = maximal_color_index;
}

void EnumerationVertexColoring::backwards() {
    determine_current_predecessors(r);
    while (!current_predecessors.empty()) {
        int i = *current_predecessors.begin();
        current_predecessors.erase(i);
        feasible_colors[ordering[i]].erase(current_solution[ordering[i]]);
        if (!feasible_colors[ordering[i]].empty()) {
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
        if (current_maximal_color < current_solution[ordering[j]]) {
            current_maximal_color = current_solution[ordering[j]];
        }
    }

    for (int j = 1; j <= current_maximal_color + 1; ++j) {
        if (j >= ub) break;
        feasible_colors_for_node.insert(j);
    }

    for (int j = 0; j < i; ++j) {
        if (graph->hasEdge(ordering[i], ordering[j])) {
            feasible_colors_for_node.erase(current_solution[ordering[j]]);
        }
    }
    if (feasible_colors.find(ordering[i]) != feasible_colors.end())
        feasible_colors.erase(ordering[i]);
    feasible_colors.insert(std::make_pair(ordering[i], feasible_colors_for_node));
}

std::vector<NetworKit::node>
    BrownsOrdinaryEnumerationVertexColoring::greedy_largest_first_ordering() {
    std::unordered_set<NetworKit::node> ordered;
    std::vector<NetworKit::node> ordering;

    NetworKit::node max_degree_node = 0;
    unsigned int max_degree = 0;
    graph->forNodes([&](NetworKit::node u) {
        if (graph->degree(u) > max_degree) {
            max_degree = graph->degree(u);
            max_degree_node = u;
        }
    });
    ordering.push_back(max_degree_node);
    ordered.insert(max_degree_node);

    std::map<NetworKit::node, std::pair<int, uint>> neighbours_in_ordering;
    int max_neighbours = 0;
    max_degree = 0;
    NetworKit::node max_neighbours_node = 0;

    graph->forNodes([&](NetworKit::node u) {
        if (u == max_degree_node) {
            return;
        }
        if (graph->hasEdge(max_degree_node, u)) {
            if (graph->degree(u) > max_degree) {
                max_degree = graph->degree(u);
                max_neighbours_node = u;
                max_neighbours = 1;
            }
            neighbours_in_ordering.insert(std::make_pair(u, std::make_pair(1, graph->degree(u))));
        } else {
            neighbours_in_ordering.insert(std::make_pair(u, std::make_pair(0, graph->degree(u))));
        }
    });

    while (ordered.size() < graph->numberOfNodes()) {
        ordering.push_back(max_neighbours_node);
        ordered.insert(max_neighbours_node);

        graph->forNodes([&](NetworKit::node v) {
            if (ordered.find(v) == ordered.end()) {
                if (graph->hasEdge(max_neighbours_node, v)) {
                    neighbours_in_ordering[v].first++;
                }
            }
        });

        max_neighbours = 0;
        max_degree = 0;
        graph->forNodes([&](NetworKit::node v) {
            if (ordered.find(v) == ordered.end()) {
                if (neighbours_in_ordering[v].first > max_neighbours) {
                    max_neighbours = neighbours_in_ordering[v].first;
                    max_degree = neighbours_in_ordering[v].second;
                    max_neighbours_node = v;
                } else if (neighbours_in_ordering[v].first == max_neighbours) {
                    if (neighbours_in_ordering[v].second > max_degree) {
                        max_degree = neighbours_in_ordering[v].second;
                        max_neighbours_node = v;
                    }
                }
            }
        });
    }
    return ordering;
}

void BrownsOrdinaryEnumerationVertexColoring::determine_current_predecessors(int r) {
    current_predecessors = std::set<int, std::greater<int>>();
    for (int i = 0; i < r; ++i) {
        current_predecessors.insert(i);
    }
}

void BrownsOrdinaryEnumerationVertexColoring::run() {
    ordering = greedy_largest_first_ordering();
    lower_bound = 1;
    upper_bound = graph->numberOfNodes();
    n = graph->numberOfNodes();

    r = 0;
    ub = upper_bound + 1;

    feasible_colors[ordering[0]].insert(1);

    while (true) {
        forwards();
        if (ub == lower_bound) {
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
    transitive_closure.resize(graph->numberOfNodes());
    for (NetworKit::node u = 0; u < graph->numberOfNodes(); ++u) {
        transitive_closure[u].resize(graph->numberOfNodes());
    }
    for (NetworKit::node u = 0; u < graph->numberOfNodes(); ++u) {
        for (NetworKit::node v = 0; v < graph->numberOfNodes(); ++v) {
            if (graph->hasEdge(u, v) && u < v) {
                transitive_closure[u][v] = true;
            }
        }
    }
    for (NetworKit::node u = 0; u < graph->numberOfNodes(); ++u) {
        for (NetworKit::node v = 0; v < graph->numberOfNodes(); ++v) {
            if (transitive_closure[u][v]) {
                for (NetworKit::node w = 0; w < graph->numberOfNodes(); ++w) {
                    if (transitive_closure[v][w]) {
                        transitive_closure[u][w] = true;
                    }
                }
            }
        }
    }
}

void ChristofidesEnumerationVertexColoring::determine_current_predecessors(int r) {
    for (int u = 0; u < n; ++u) {
        if (transitive_closure[u][ordering[r]]) {
            current_predecessors.insert(u);
        }
    }
}

const std::vector<std::vector<bool>>& ChristofidesEnumerationVertexColoring::getTransitiveClosure()
const {
    assureFinished();
    return transitive_closure;
}

void ChristofidesEnumerationVertexColoring::run() {
    calculate_transitive_closure();
    graph->forNodes([&](NetworKit::node u) {
        ordering.push_back(u);
    });
    lower_bound = 1;
    upper_bound = graph->numberOfNodes();
    n = graph->numberOfNodes();

    r = 0;
    ub = upper_bound + 1;

    feasible_colors[ordering[0]].insert(1);

    while (true) {
        forwards();
        if (ub == lower_bound) {
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
    std::vector<NetworKit::node>& subgraph,
    NetworKit::node new_node
) {
    std::unordered_set<NetworKit::node> visited;
    int number_of_neighbours_in_component = 0;
    for (NetworKit::node u : subgraph) {
        if (visited.find(u) == visited.end()) {
            std::vector<NetworKit::node> component;
            std::queue<NetworKit::node> queue;
            number_of_neighbours_in_component = 0;
            queue.push(u);
            visited.insert(u);
            while (!queue.empty()) {
                NetworKit::node v = queue.front();
                queue.pop();
                if (graph->hasEdge(new_node, v)) {
                    ++number_of_neighbours_in_component;
                }
                component.push_back(v);
                for (NetworKit::node w : subgraph) {
                    if (graph->hasEdge(v, w) && visited.find(w) == visited.end()) {
                        queue.push(w);
                        visited.insert(w);
                    }
                }
            }
            if (number_of_neighbours_in_component == 1) {
                return component;
            }
        }
    }
    return std::vector<NetworKit::node>();
}

bool BrelazEnumerationVertexColoring::is_interchangeable(
    std::vector<int>& K, NetworKit::node new_node
) {
    for (int alpha : K) {
        for (int beta : K) {
            if (alpha == beta) break;
            std::vector<NetworKit::node> subgraph;
            for (const auto& [node, color] : current_solution) {
                if (color == alpha || color == beta) {
                    subgraph.push_back(node);
                }
            }
            auto interchangeable_component = interchange_component(subgraph, new_node);
            if (interchangeable_component.size() > 0) {
                for (NetworKit::node v : interchangeable_component) {
                    current_solution[v] = (alpha == current_solution[v]) ? beta : alpha;
                    if (graph->hasEdge(new_node, v)) {
                        current_solution[new_node] = (alpha == current_solution[v]) ? beta : alpha;
                    }
                }
                return true;
            }
        }
    }
    return false;
}

std::vector<NetworKit::node>
BrelazEnumerationVertexColoring::saturation_largest_first_with_interchange() {
    std::vector<NetworKit::node> ordering;
    NetworKit::node max_degree_node = 0;
    uint max_degree = 0;

    graph->forNodes([&](NetworKit::node u) {
        if (graph->degree(u) > max_degree) {
            max_degree = graph->degree(u);
            max_degree_node = u;
        }
    });

    current_solution.insert(std::make_pair(max_degree_node, 1));
    ordering.push_back(max_degree_node);

    auto satur_comp = [&](
        std::tuple<int, int, NetworKit::node> a, std::tuple<int, int, NetworKit::node> b) {
            if (std::get<0>(a) == std::get<0>(b)) {
            if (std::get<1>(a) == std::get<1>(b)) {
                return std::get<2>(a) < std::get<2>(b);
            }
            return std::get<1>(a) > std::get<1>(b);
        }
        return std::get<0>(a) > std::get<0>(b);
    };

    auto saturation =
        std::set <std::tuple<int, int, NetworKit::node>, decltype(satur_comp)>(satur_comp);

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

    int max_color = 1;
    int max_clique_size = 0;

    while (current_solution.size() < n) {
        auto max_saturation = saturation.begin();
        NetworKit::node u = std::get<2>(*max_saturation);
        saturation.erase(max_saturation);
        ordering.push_back(u);

        int first_valid_color = 1;
        std::set<int> forbidden_colors;
        std::vector<std::pair<int, std::vector<NetworKit::node>>> neighbours_color_count;
        neighbours_color_count.resize(max_color + 1);
        graph->forNeighborsOf(u, [&](NetworKit::node v) {
            if (current_solution.find(v) != current_solution.end()) {
                forbidden_colors.insert(current_solution[v]);
                if (neighbours_color_count[current_solution[v]].first == 0) {
                    neighbours_color_count[current_solution[v]].first = 1;
                    neighbours_color_count[current_solution[v]].second.push_back(v);
                }
            }
        });
        while (forbidden_colors.find(first_valid_color) != forbidden_colors.end()) {
            ++first_valid_color;
        }

        if (first_valid_color <= max_color) {
            current_solution.insert(std::make_pair(u, first_valid_color));
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
            auto interchange = is_interchangeable(K, u);
            if (!interchange) {
                current_solution.insert(std::make_pair(u, ++max_color));
            } else {
                if (max_clique_size == 0) {
                    max_clique_size = max_color;
                }
            }
        }
        for (const auto & [saturation_degree, degree, node] : saturation) {
            if (graph->hasEdge(u, node)) {
                neighbours_colors[node].insert(current_solution[u]);
                saturation.erase(std::make_tuple(saturation_degree, degree, node));
                saturation.insert(
                    std::make_tuple(neighbours_colors[node].size(), degree - 1, node));
            }
        }
    }

    upper_bound = max_color;
    lower_bound = max_clique_size;

    best_solution = current_solution;

    return ordering;
}

std::vector<int> BrelazEnumerationVertexColoring::get_representatives_of_adjacent_predecessors(
    int i
) {
    std::vector<int> representatives(ub, INT_MAX);
    for (int j = 0; j < i; ++j) {
        if (graph->hasEdge(ordering[j], ordering[i]) && current_solution[ordering[j]] < ub) {
            if (representatives[current_solution[ordering[j]]] > j) {
                representatives[current_solution[ordering[j]]] = j;
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
        feasible_colors[ordering[i]].erase(current_solution[ordering[i]]);
        if (!feasible_colors[ordering[i]].empty()) {
            r = i;
            return;
        }
    }
    r = 0;
}

void BrelazEnumerationVertexColoring::run() {
    n = graph->numberOfNodes();

    ordering = saturation_largest_first_with_interchange();

    if (lower_bound != upper_bound) {
        r = 0;
        ub = upper_bound;

        feasible_colors[ordering[0]].insert(1);

        while (true) {
            forwards();
            if (ub == lower_bound) {
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
} /* namespace Koala */
