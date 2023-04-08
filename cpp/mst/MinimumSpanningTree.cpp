/*
 * MinimumSpanningTree.cpp
 *
 *  Created on: 07.04.2023
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include <mst/MinimumSpanningTree.hpp>

#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/structures/UnionFind.hpp>

#include "mst/BoruvkaMST.hpp"
#include "mst/KktMST.hpp"

namespace Koala {

MinimumSpanningTree::MinimumSpanningTree(
    NetworKit::Graph &graph) : graph(std::make_optional(graph)) { }

const NetworKit::Graph& MinimumSpanningTree::getForest() const {
    assureFinished();
    return *tree;
}

void KruskalMinimumSpanningTree::run() {
    hasRun = true;
    std::vector<NetworKit::WeightedEdge> sorted_edges(
        graph->edgeWeightRange().begin(), graph->edgeWeightRange().end());
    Aux::Parallel::sort(sorted_edges.begin(), sorted_edges.end());
    tree = std::make_optional(NetworKit::GraphTools::copyNodes(*graph));
    NetworKit::UnionFind union_find(graph->upperNodeIdBound());
    for (const auto &e : sorted_edges) {
        if (union_find.find(e.u) != union_find.find(e.v)) {
            tree->addEdge(e.u, e.v, e.weight);
            union_find.merge(e.u, e.v);
        }
    }
}

void BoruvkaMinimumSpanningTree::run() {
    hasRun = true;
    auto algorithm = MST::BoruvkaMST(*graph);
    tree = std::make_optional(algorithm.getSpanningTree());
}

void KargerKleinTarjanMinimumSpanningTree::run() {
    hasRun = true;
    auto algorithm = MST::KktMST(*graph);
    tree = std::make_optional(algorithm.getSpanningTree());
}

} /* namespace Koala */
