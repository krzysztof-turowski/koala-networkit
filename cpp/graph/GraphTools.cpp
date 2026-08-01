/*
 * GraphTools.cpp
 *
 *  Created on: 12.11.2021
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include <graph/GraphTools.hpp>

#include <algorithm>
#include <iterator>
#include <list>
#include <queue>
#include <random>
#include <set>
#include <stdexcept>
#include <unordered_set>
#include <vector>

#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/graph/AdjListGraph.hpp>

namespace Koala {

namespace GraphTools {

NetworKit::Graph toComplement(const NetworKit::Graph &G) {
    NetworKit::Graph GC(G.upperNodeIdBound(), false, false);
    for (NetworKit::node v = 0; v < G.upperNodeIdBound(); v++) {
        if (G.hasNode(v)) {
            std::set<NetworKit::node> neighbors(
                G.neighborRange(v).begin(), G.neighborRange(v).end());
            GC.forNodes([&](NetworKit::node u) {
                if (u < v && !neighbors.count(u)) {
                    GC.addEdge(u, v);
                }
            });
        } else {
            GC.removeNode(v);
        }
    }
    return GC;
}

NetworKit::Graph convertDirectedGraphToUndirected(const NetworKit::Graph &G, bool weighted) {
    NetworKit::Graph Gout(G.upperNodeIdBound(), weighted, false);
    for (NetworKit::node u = 0; u < G.upperNodeIdBound(); ++u) {
        if (!G.hasNode(u)) {
            Gout.removeNode(u);
        }
    }
    G.forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) {
        Gout.increaseWeight(u, v, w);
    });
    Gout.removeMultiEdges();
    return Gout;
}

NetworKit::Graph convertUndirectedGraphToDirected(const NetworKit::Graph &G, bool weighted) {
    NetworKit::Graph Gout(G.upperNodeIdBound(), weighted, true);
    for (NetworKit::node u = 0; u < G.upperNodeIdBound(); ++u) {
        if (!G.hasNode(u)) {
            Gout.removeNode(u);
        }
    }
    G.forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) {
        Gout.increaseWeight(u, v, w);
        Gout.increaseWeight(v, u, w);
    });
    Gout.removeMultiEdges();
    return Gout;
}

void assureUndirectedGraph(const NetworKit::Graph &G) {
    if (G.isDirected()) {
        throw std::runtime_error("Algorithm requires an undirected graph.");
    }
}

void ensureDirectedGraph(NetworKit::Graph &G) {
    if (!G.isDirected()) {
        G = NetworKit::Graph(G, G.isWeighted(), true);
    }
}

bool isConnected(const NetworKit::Graph &G) {
    NetworKit::Graph Gout(G);
    auto connected_components = NetworKit::ConnectedComponents(Gout);
    connected_components.run();
    return connected_components.getComponents().size() == 1;
}

std::vector<std::vector<NetworKit::node>> complementComponents(
        const NetworKit::Graph &G, const std::vector<NetworKit::node> &nodes) {
    std::list<NetworKit::node> unvisited(nodes.begin(), nodes.end());
    std::vector<bool> in_unvisited(G.upperNodeIdBound());
    std::vector<std::list<NetworKit::node>::iterator> position(G.upperNodeIdBound());
    for (auto it = unvisited.begin(); it != unvisited.end(); ++it) {
        in_unvisited[*it] = true;
        position[*it] = it;
    }
    std::vector<std::vector<NetworKit::node>> components;
    std::queue<NetworKit::node> queue;
    std::vector<NetworKit::node> neighbours;

    while (!unvisited.empty()) {
        components.emplace_back();
        queue.push(unvisited.front());
        in_unvisited[unvisited.front()] = false;
        unvisited.pop_front();
        while (!queue.empty()) {
            auto u = queue.front();
            queue.pop();
            components.back().push_back(u);

            neighbours.clear();
            for (auto v : G.neighborRange(u)) {
                if (v < in_unvisited.size() && in_unvisited[v]) {
                    neighbours.push_back(v);
                    in_unvisited[v] = false;
                    unvisited.erase(position[v]);
                }
            }

            while (!unvisited.empty()) {
                auto v = unvisited.front();
                queue.push(v);
                in_unvisited[v] = false;
                unvisited.pop_front();
            }

            for (auto v : neighbours) {
                in_unvisited[v] = true;
                unvisited.push_back(v);
                position[v] = std::prev(unvisited.end());
            }
        }
    }

    return components;
}

bool hasMultiEdges(const NetworKit::Graph &G) {
    std::unordered_set<NetworKit::Edge> edges;
    for (const auto &edge : G.edgeRange()) {
        if (!edges.insert(edge).second) {
            return true;
        }
    }
    return false;
}

NetworKit::Graph makeConnected(const NetworKit::Graph &G) {
    NetworKit::Graph Gout(G);
    auto connected_components = NetworKit::ConnectedComponents(Gout);
    connected_components.run();
    auto components = connected_components.getComponents();
    for (NetworKit::count i = 1; i < connected_components.numberOfComponents(); i++) {
        Gout.addEdge(components[0][0], components[i][0]);
    }
    return Gout;
}

bool hasDistinctIntegerWeights(const NetworKit::Graph &G) {
    std::set<NetworKit::edgeweight> weights;
    G.forEdges([&](NetworKit::node, NetworKit::node, NetworKit::edgeweight w) {
        if (w == static_cast<NetworKit::count>(w)) {
            weights.insert(w);
        }
    });
    return weights.size() == G.numberOfEdges();
}

NetworKit::Graph makeDistinctIntegerWeights(
        const NetworKit::Graph &G, NetworKit::edgeweight max_weight) {
    if (G.numberOfEdges() > max_weight) {
        throw std::runtime_error(
            "Operation impossible: more edges than available distinct weights");
    }

    NetworKit::Graph Gout(G, true, G.isDirected(), G.hasEdgeIds());
    if (G.numberOfEdges() == 0) {
        return Gout;
    }

    auto max_weight_integer = static_cast<NetworKit::count>(max_weight);
    std::random_device device;
    std::mt19937_64 generator(device());
    std::unordered_set<NetworKit::count> selected_weights;
    selected_weights.reserve(G.numberOfEdges());

    for (NetworKit::count value = max_weight_integer - G.numberOfEdges() + 1;
            value <= max_weight_integer; ++value) {
        std::uniform_int_distribution<NetworKit::count> distribution(1, value);
        auto candidate = distribution(generator);
        if (selected_weights.contains(candidate)) {
            selected_weights.insert(value);
        } else {
            selected_weights.insert(candidate);
        }
    }

    std::vector<NetworKit::count> weights(selected_weights.begin(), selected_weights.end());
    std::shuffle(weights.begin(), weights.end(), generator);
    NetworKit::count id = 0;
    Gout.forEdges([&](NetworKit::node u, NetworKit::node v) {
        Gout.setWeight(u, v, weights[id++]);
    });
    return Gout;
}

}  // namespace GraphTools

}  // namespace Koala
