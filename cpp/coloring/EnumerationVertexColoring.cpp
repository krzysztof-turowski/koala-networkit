/*
* EnumerationVertexColoring.cpp
*
* Created on: 06.02.2023
*   Author: Zofia Glapa (zofia.glapa@student.uj.edu.pl)
*/

#include <coloring/EnumerationVertexColoring.hpp>

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
    for (int i = 0; i < r; ++i) {
        current_predecessors.push(i);
    }
    while (!current_predecessors.empty()) {
        int i = current_predecessors.top();
        current_predecessors.pop();
        feasible_colors[ordering[i]].erase(current_solution[ordering[i]]);
        if (!feasible_colors[ordering[i]].empty()) {
            r = i;
            return;
        }
    }
    r = 0;
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

void BrownsOrdinaryEnumerationVertexColoring::determine_feasible_colors(int i) {
    std::set<int> feasible_colors_for_node;
    int current_maximal_color = 0;
    for (int j = 0; j < i; ++j) {
        if (current_maximal_color < current_solution[ordering[j]]) {
            current_maximal_color = current_solution[ordering[j]];
        }
    }

    for (int j = 1; j <= current_maximal_color + 1; ++j) {
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

} /* namespace Koala */
