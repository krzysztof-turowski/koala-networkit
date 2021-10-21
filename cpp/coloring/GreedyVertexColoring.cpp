/*
 * GreedyVertexColoring.cpp
 *
 *  Created on: 22.10.2021
 *      Author: Krzysztof Turowski
 */

#include <coloring/GreedyVertexColoring.hpp>

namespace Koala {

int GreedyVertexColoring::randomSequential(
        const NetworKit::Graph &G, std::map<NetworKit::node, int> &colors) {
    int max_color = 0;
    G.forNodesInRandomOrder([&](NetworKit::node v) {
        int color = greedy_color(G, v, colors)->second;
        if (color > max_color) {
            max_color = color;
        }
    });
    return max_color;
}

int GreedyVertexColoring::largestFirst(
        const NetworKit::Graph &G, std::map<NetworKit::node, int> &colors) {
    std::vector<NetworKit::node> vertices(largest_first_ordering(G, colors));
    int max_color = 0;
    for (NetworKit::node v : vertices) {
        int color = greedy_color(G, v, colors)->second;
        if (color > max_color) {
            max_color = color;
        }
    }
    return max_color;
}

std::map<NetworKit::node, int>::iterator GreedyVertexColoring::greedy_color(
        const NetworKit::Graph &G, NetworKit::node v, std::map<NetworKit::node, int> &colors) {
    std::map<NetworKit::node, int>::iterator first(colors.find(v));
    if (first != colors.end()) {
        return first;
    }

    int max_color = G.degree(v) + 2;
    std::vector<bool> forbidden(max_color, false);
    forbidden[0] = true;
    for (const NetworKit::node u : G.neighborRange(v)) {
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

std::vector<NetworKit::node> GreedyVertexColoring::largest_first_ordering(
        const NetworKit::Graph &G, std::map<NetworKit::node, int> &colors) {
    std::vector<NetworKit::node> vertices;
    G.forNodes([&](NetworKit::node v) {
        if (colors.find(v) == colors.end()) {
            vertices.push_back(v);
        }
    });
    std::sort(
        vertices.begin(), vertices.end(),
        [&](auto &u, auto &v) { return G.degree(u) > G.degree(v);
    });
    return vertices;
}

} /* namespace Koala */
