#include "KktMST.hpp"
#include "MinimalSpanningTree.hpp"
#include "BoruvkaMST.hpp"
#include "MSTV.hpp"
#include "utils.hpp"
#include <iostream>
#include <networkit/components/ConnectedComponents.hpp>

namespace MST {

    Graph discardEdgesRandomlyButKeepConnected(const Graph& G) {
        /// returned graph is still always connected - by edges of weight +inf, if needed.
        srand(time(NULL));
        Graph H(G, true, false);
        G.forEdges([&](const node u, const node v) {
           if (rand() & 1) {
                H.removeEdge(u, v);
           }
        });
        makeConnected(H, std::numeric_limits<NetworKit::edgeweight>::max()); // doesn't invalidate runtime analysis!
        return H;
    }

    KktMST::KktMST(const Graph& G, bool storeOriginalGraph) : MinimalSpanningTree(G, storeOriginalGraph) {
        if (G.numberOfEdges() == 0) {
            return;
        }
        auto G_prim_MST = BoruvkaMST(G, false, false, true);
        G_prim_MST.getSpanningTree().forEdges([&](const node u, const node v, edgeweight ew) {
           mst.addEdge(u, v, ew);
        });
        if (G_prim_MST.contractedG_iAfterLimitedSteps.value().numberOfEdges() == 0)
            return;

        Graph H = discardEdgesRandomlyButKeepConnected(G_prim_MST.contractedG_iAfterLimitedSteps.value());
        KktMST H_kktMSt(H, false);
        const Graph& F = H_kktMSt.getSpanningTree();
        removeF_HeavyEdges(G_prim_MST.contractedG_iAfterLimitedSteps.value(), F);

        auto recursiveCallOutput = KktMST(G_prim_MST.contractedG_iAfterLimitedSteps.value());
        recursiveCallOutput.getSpanningTree().forEdges([&](const node u, const node v, const edgeweight ew) {
            auto eidInContractedG_i = G_prim_MST.contractedG_iAfterLimitedSteps.value().edgeId(u, v);
            auto[origU, origV] = G_prim_MST.contractedG_iGetOriginalNodesFromEdgeId.value()[eidInContractedG_i];
            mst.addEdge(origU, origV, ew);
        });
    }
}
