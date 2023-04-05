#include "KruskalMST.hpp"
#include "MinimalSpanningTree.hpp"
#include "utils.hpp"
#include <networkit/graph/Graph.hpp>
#include <networkit/structures/UnionFind.hpp>

namespace MST {

    KruskalMST::KruskalMST(const Graph& G, bool storeOriginalGraph) : MinimalSpanningTree(G, storeOriginalGraph) {
        auto sets = NetworKit::UnionFind(G.numberOfNodes());
        auto edgeRange = G.edgeWeightRange();
        std::vector<NetworKit::WeightedEdge> edges(edgeRange.begin(), edgeRange.end());
        std::sort(edges.begin(), edges.end(), [](const auto& lhs, const auto& rhs) { return lhs.weight < rhs.weight; });
        for (auto& e: edges) {
            if (sets.find(e.u) != sets.find(e.v)) {
                mst.addEdge(e.u, e.v, e.weight);
                sets.merge(e.u, e.v);
            }
        }
    }
}
