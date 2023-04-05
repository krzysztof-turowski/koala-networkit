#pragma once

#include <networkit/graph/Graph.hpp>
#include <networkit/components/ConnectedComponents.hpp>

#define MST_DEBUG

namespace MST {
    using NetworKit::Graph;
    using NetworKit::node;
    using NetworKit::edgeid;
    using NetworKit::index;
    using NetworKit::count;
    using NetworKit::edgeweight;
    using NodePair = std::pair<node, node>;

    inline constexpr count log2(count n) {
        // returns floor of log2(n).
        return n > 1 ? 1 + log2(n >> 1) : 0;
    }

    inline node theOtherEnd(const NetworKit::Edge& e, node u) {
        return e.u == u ? e.v : e.u;
    }

    inline std::ostream& operator<<(std::ostream& os, const Graph& G) {
        os << "Graph Info. N: " << G.numberOfNodes() << " M: " << G.numberOfEdges()  << std::endl;
        G.forEdges([&](const node u, const node v, edgeweight ew, edgeid eid) {
            os << u << " -> " << v << ", w: " << ew << ", eid: " << eid << std::endl;
        });
        os << "Total weight: " << G.totalEdgeWeight() << std::endl;
        return os;
    }

    inline void makeConnected(Graph& G, const edgeweight ew) {
        auto ccs = NetworKit::ConnectedComponents(G);
        ccs.run();
        auto components = ccs.getComponents();
        for(index i = 1; i < ccs.numberOfComponents(); i++) {
            G.addEdge(components[0][0], components[i][0], ew);
        }
    }

    namespace WeightedUndirectedEdgeDecoupled {
        /**
         * `type.first = node` of an edge (only one, because that's what Graph.weightNeighborRange() exposes).
         * **/
        using type = std::pair<node, edgeweight>;

        inline bool less(const type& a, const type& b) {
            return a.second < b.second;
        }
    }
}
