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

NetworKit::Graph convertDirectedGraphToUndirected(NetworKit::Graph &G, bool weighted) {
    NetworKit::Graph Gout(G.numberOfNodes(), weighted, false);
    G.forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) {
        Gout.addEdge(u, v, w);
    });
    Gout.removeMultiEdges();
    return Gout;
}

NetworKit::Graph convertUndirectedGraphToDirected(NetworKit::Graph &G, bool weighted) {
    NetworKit::Graph Gout(G.numberOfNodes(), weighted, true);
    G.forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) {
        Gout.addEdge(u, v, w);
        Gout.addEdge(v, u, w);
    });
    Gout.removeMultiEdges();
    return Gout;
}

}  // namespace GraphTools

}  // namespace Koala
