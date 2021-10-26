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
        [&](auto &u, auto &v) { return graph->degree(u) > graph->degree(v);
    });
    return vertices;
}

void SmallestLastVertexColoring::run() {
    std::vector<NetworKit::node> vertices(smallest_last_ordering());
    int max_color = 0;
    for (const auto v : vertices) {
        int color = greedy_color(v)->second;
        if (color > max_color) {
            max_color = color;
        }
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

} /* namespace Koala */
