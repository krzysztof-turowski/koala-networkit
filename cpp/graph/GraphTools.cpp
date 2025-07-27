/*
 * GraphTools.cpp
 *
 *  Created on: 12.11.2021
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include <set>

#include <networkit/graph/Graph.hpp>

#include <graph/GraphTools.hpp>

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

void printGraph(const NetworKit::Graph& graph) {
    std::cout << "Graph structure" << std::endl;
    for (auto edge : graph.edgeWeightRange()) {
        std::cout << edge.u << " - " << edge.v << " weight: " << edge.weight << std::endl;
    }
}

NetworKit::Graph convertDirectedGraphToUndirected(NetworKit::Graph& graph) {
    NetworKit::Graph result(graph.numberOfNodes());
    for (auto [u, v] : graph.edgeRange()) {
        result.addEdge(u, v);
    }
    result.removeMultiEdges();
    return result;
}

}  // namespace GraphTools

}  // namespace Koala
