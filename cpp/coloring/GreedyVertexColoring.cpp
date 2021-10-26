/*
 * GreedyVertexColoring.cpp
 *
 *  Created on: 22.10.2021
 *      Author: Krzysztof Turowski
 */

#include <coloring/GreedyVertexColoring.hpp>

namespace Koala {

GreedyVertexColoring::GreedyVertexColoring(
        const NetworKit::Graph &graph) : graph(std::make_optional(graph)) { }

const std::map<NetworKit::node, int>& GreedyVertexColoring::getColoring() const {
    assureFinished();
    return colors;
}

std::map<NetworKit::node, int>::iterator GreedyVertexColoring::greedy_color(NetworKit::node v) {
    std::map<NetworKit::node, int>::iterator first(colors.find(v));
    if (first != colors.end()) {
        return first;
    }

    int max_color = graph->degree(v) + 2;
    std::vector<bool> forbidden(max_color, false);
    forbidden[0] = true;
    for (const auto u : graph->neighborRange(v)) {
        std::map<NetworKit::node, int>::iterator second(colors.find(u));
        if (second != colors.end()) {
            int color = second->second;
            if (color < max_color) {
                forbidden[color] = true;
            }
        }
    }
    int color = std::find(forbidden.begin(), forbidden.end(), false) - forbidden.begin();
    return colors.insert(std::make_pair(v, color)).first;
}

void RandomSequentialVertexColoring::run() {
    int max_color = 0;
    graph->forNodesInRandomOrder([&](NetworKit::node v) {
        int color = greedy_color(v)->second;
        if (color > max_color) {
            max_color = color;
        }
    });
    hasRun = true;
}

void LargestFirstVertexColoring::run() {
    std::vector<NetworKit::node> vertices(largest_first_ordering());
    int max_color = 0;
    for (const auto v : vertices) {
        int color = greedy_color(v)->second;
        if (color > max_color) {
            max_color = color;
        }
    }
    hasRun = true;
}

std::vector<NetworKit::node> LargestFirstVertexColoring::largest_first_ordering() {
    std::vector<NetworKit::node> vertices;
    graph->forNodes([&](NetworKit::node v) {
        if (colors.find(v) == colors.end()) {
            vertices.push_back(v);
        }
    });
    std::sort(
        vertices.begin(), vertices.end(),
        [&](auto &u, auto &v) { return graph->degree(u) > graph->degree(v); }
    );
    return vertices;
}

} /* namespace Koala */
