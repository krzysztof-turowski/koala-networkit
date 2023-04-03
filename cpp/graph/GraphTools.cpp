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
    NetworKit::Graph GC(G.numberOfNodes(), false, false);
    GC.forNodes([&](NetworKit::node v) {
        if (G.hasNode(v)) {
            std::set<NetworKit::node> neighbors(G.neighborRange(v).begin(), G.neighborRange(v).end());
            GC.forNodes([&](NetworKit::node u) {
                if (G.hasNode(u) && u < v && !neighbors.count(u)) {
                    GC.addEdge(u, v);
                }
            });
        }
    });
    return GC;
}

}  // namespace GraphTools

}  // namespace Koala
