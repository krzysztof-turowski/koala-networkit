/*
 * GreedyVertexColoring.cpp
 *
 *  Created on: 22.10.2021
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include <coloring/GreedyVertexColoring.hpp>

#include <networkit/auxiliary/BucketPQ.hpp>
#include <networkit/graph/GraphTools.hpp>

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
    graph->forInNeighborsOf(v, [&](NetworKit::node u) {
        std::map<NetworKit::node, int>::iterator second(colors.find(u));
        if (second != colors.end()) {
            int color = second->second;
            if (color < max_color) {
                forbidden[color] = true;
            }
        }
    });
    int color = std::find(forbidden.begin(), forbidden.end(), false) - forbidden.begin();
    return colors.insert(std::make_pair(v, color)).first;
}

void RandomSequentialVertexColoring::run() {
    graph->forNodesInRandomOrder([&](NetworKit::node v) {
        greedy_color(v)->second;
    });
    hasRun = true;
}

void LargestFirstVertexColoring::run() {
    std::vector<NetworKit::node> vertices(largest_first_ordering());
    for (const auto v : vertices) {
        greedy_color(v)->second;
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
        [&](auto &u, auto &v) { return graph->degree(u) > graph->degree(v);
    });
    return vertices;
}

void SmallestLastVertexColoring::run() {
    std::vector<NetworKit::node> vertices(smallest_last_ordering());
    for (const auto v : vertices) {
        greedy_color(v)->second;
    }
    hasRun = true;
}

std::vector<NetworKit::node> SmallestLastVertexColoring::smallest_last_ordering() {
    NetworKit::count max_degree = NetworKit::GraphTools::maxDegree(*graph);
    Aux::BucketPQ queue(graph->numberOfNodes(), 0, max_degree);
    graph->forNodes([&](NetworKit::node v) {
        if (colors.find(v) == colors.end()) {
            queue.insert(graph->degree(v), v);
        }
    });

    std::vector<NetworKit::node> vertices(queue.size());
    while (!queue.empty()) {
        auto v = queue.extractMin().second;
        vertices[queue.size()] = v;
        graph->forInNeighborsOf(v, [&](NetworKit::node u) {
            if (queue.contains(u)) {
                queue.changeKey(queue.getKey(u) - 1, u);
            }
        });
    }
    return vertices;
}

void SaturatedLargestFirstVertexColoring::run() {
    NetworKit::count max_degree = NetworKit::GraphTools::maxDegree(*graph);
    std::map<NetworKit::node, std::set<int>> saturations;
    graph->forEdges([&](NetworKit::node u, NetworKit::node v) {
        if (colors.find(u) == colors.end() && colors.find(v) != colors.end()) {
            saturations[u].insert(colors[v]);
        } else if (colors.find(u) != colors.end() && colors.find(v) == colors.end()) {
            saturations[v].insert(colors[u]);
        }
    });

    NetworKit::count limit = (max_degree + 1) * max_degree;
    Aux::BucketPQ queue(graph->numberOfNodes(), -limit, 0);
    graph->forNodes([&](NetworKit::node v) {
        if (colors.find(v) == colors.end()) {
            queue.insert(-saturations[v].size() * max_degree - graph->degree(v), v);
        }
    });

    while (!queue.empty()) {
        auto v = queue.extractMin().second;
        int color = greedy_color(v)->second;
        graph->forInNeighborsOf(v, [&](NetworKit::node u) {
            if (queue.contains(u) && !saturations[u].count(color)) {
                saturations[u].insert(color);
                queue.changeKey(queue.getKey(u) - max_degree, u);
            }
        });
    }
    hasRun = true;
}

} /* namespace Koala */
