#include "BoruvkaMST.hpp"
#include "MinimalSpanningTree.hpp"
#include "utils.hpp"
#include <networkit/graph/Graph.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <memory>
#include <optional>

namespace MST {

    // Note that this constructor merely begins to initialize FullBranchingTree. It will be fully constructed
    // only after BoruvkaMST() returns.
    FullBranchingTree::FullBranchingTree(const Graph& G) : leaves(G.upperNodeIdBound()),
                                                           tree(G.upperNodeIdBound(), true, true) {}

    node FullBranchingTree::getRoot() const { return tree.upperNodeIdBound() - 1; }

    const Graph& FullBranchingTree::getTree() const { return tree; }

    std::optional<node> FullBranchingTree::getParent(node u) const {
        if (u == getRoot()) {
            return {std::nullopt};
        }
        // make use of the fact that there is only one incoming edge.
        getTree().forInNeighborsOf(u, [&u](NetworKit::node, NetworKit::node parent, NetworKit::edgeweight) {
            u = parent;
        });
        return {u};
    }

    class Bookkeeping {
        /** Helper for Boruvka algorithm + responsible for populating fullBranchingTree. **/
    private:
        std::vector<node> V_iToV_i_plus_1;
        count G_i_plus_1_nodeCount = 0;
        Graph F_i;
        // fullBranchingTree structures
        std::optional<FullBranchingTree>& fullBranchingTree;
        std::vector<node> V_iToFullBranchingTreeNode;


    public:
        std::vector<NodePair> E_iToE_0;

        explicit Bookkeeping(const Graph& G, std::optional<FullBranchingTree>& fullBranchingTree)
                : E_iToE_0(G.upperEdgeIdBound()), V_iToV_i_plus_1(G.upperNodeIdBound()),
                  F_i(G.upperNodeIdBound(), true, false),
                  fullBranchingTree(fullBranchingTree) {
            /** we assume that G.has indexed edges. **/
//        G.parallelForEdges([&](node u, node v, edgeid edgeId) {
            G.forEdges([&](node u, node v, edgeweight, edgeid edgeId) {
                E_iToE_0[edgeId] = {u, v};
            });
            if (fullBranchingTree) {
                V_iToFullBranchingTreeNode.resize(G.upperNodeIdBound());
                std::iota(V_iToFullBranchingTreeNode.begin(), V_iToFullBranchingTreeNode.end(), 0);
            }
        }

        [[nodiscard]]
        NodePair getOriginalNodes(edgeid edgeId) const {
            return E_iToE_0.at(edgeId);
        }


        void addEdgeToF_i(node u, node v, edgeweight edgeWeight) {
            if (!F_i.hasEdge(u, v)) {
                F_i.addEdge(u, v, edgeWeight);
            }
        }

        //TODO: consider parallelism here.
        void updateFullBranchingTree(const std::vector<std::vector<node>>& connectedComponents,
                                     const std::vector<WeightedUndirectedEdgeDecoupled::type>& edgesForMst) {
            std::vector<node> V_i_plus_1ToFullBranchingTreeNode;
            for (const auto& cc: connectedComponents) {
                auto parentNode = fullBranchingTree->tree.addNode();
                for (auto u: cc) {
                    fullBranchingTree->tree.addEdge(parentNode, V_iToFullBranchingTreeNode[u], edgesForMst[u].second);
                }
                V_i_plus_1ToFullBranchingTreeNode.push_back(parentNode);
            }
            V_iToFullBranchingTreeNode = std::move(V_i_plus_1ToFullBranchingTreeNode);
        }

        void contractF_iNodes(const std::vector<std::vector<node>>& connectedComponents) {
            /** Traverses F_i to obtain new indices. **/
            for (count i = 0; i < connectedComponents.size(); i++) {
                for (auto u: connectedComponents[i]) {
                    V_iToV_i_plus_1[u] = i;
                }
            }
        }

        Graph contractG_iTo_G_i_plus_1(const Graph& G_i, const std::vector<std::vector<node>>& connectedComponents) {
            /** Eisner: Chapter 4.12 **/
            contractF_iNodes(connectedComponents);
            G_i_plus_1_nodeCount = connectedComponents.size();
            std::vector<std::vector<NetworKit::WeightedEdgeWithId>> bins1(G_i_plus_1_nodeCount);
            G_i.forEdges([&](node a, node b, edgeweight w, edgeid edgeId) {
                node v = std::max(V_iToV_i_plus_1[a], V_iToV_i_plus_1[b]);
                node u = std::min(V_iToV_i_plus_1[a], V_iToV_i_plus_1[b]);
                bins1[v].emplace_back(u, v, w, edgeId);
            });
            // we won't use G_i anymore

            std::vector<std::vector<NetworKit::WeightedEdgeWithId>> bins2(G_i_plus_1_nodeCount);
            for (node v = G_i_plus_1_nodeCount; v-- > 0;) { // iterates through [G_i_plus_1_nodeCount, 0]
                for (const auto& e: bins1[v]) {
                    node u = theOtherEnd(e, v);
                    bins2[u].emplace_back(e);
                }
            }
            /**
             * bins2[u] now contains edges sorted in lexicographic order by their endpoints.
             * Hence, edges with the same endpoints are stored consecutively.
             */
            Graph G_i_plus_1(G_i_plus_1_nodeCount, true, false);
            G_i_plus_1.indexEdges(true);
            std::vector<NodePair> E_i_plus_1ToE_0;
            for (node u = 0; u < G_i_plus_1_nodeCount; ++u) {
                for (auto left = bins2[u].begin(), right = bins2[u].begin(); left != bins2[u].end(); left = right) {
                    node endPoint = theOtherEnd(*left, u);
                    while ((right != bins2[u].end()) && (theOtherEnd(*right, u) == endPoint)) {
                        right = std::next(right);
                    }
                    /* [Left, Right) is a range of edges with the same endpoint. */
                    /* if (not a self loop) && (`endpoint` hasn't already added  this edge) */
                    if (endPoint != u && !G_i_plus_1.hasEdge(u, endPoint)) {
                        auto lightestEdge = std::min_element(left, right, [](const auto& a, const auto& b) {
                            return a.weight < b.weight;
                        });
                        G_i_plus_1.addEdge(u, endPoint, lightestEdge->weight);
                        E_i_plus_1ToE_0.emplace_back(getOriginalNodes(lightestEdge->eid));
                    }
                }
            }
            // work is done: prepare for the next iteration.
            F_i = Graph(G_i_plus_1_nodeCount, true, false);
            E_iToE_0 = std::move(E_i_plus_1ToE_0);
            return G_i_plus_1;
        }


        Graph prepareForNextIteration(const Graph& G_i,
                                      const std::vector<WeightedUndirectedEdgeDecoupled::type>& edgesForMst) {
            auto ccs = NetworKit::ConnectedComponents(F_i);
            ccs.run();
            auto connectedComponents = ccs.getComponents();
            if (fullBranchingTree) {
                updateFullBranchingTree(connectedComponents, edgesForMst);
            }
            return contractG_iTo_G_i_plus_1(G_i, connectedComponents);
        }
    };


    /** Constructor initializes `this->mst` and potentially `this->fullBranchingTree`. **/
    BoruvkaMST::BoruvkaMST(const Graph& G, bool storeOriginalGraph, bool storeFullBranchingTree, bool storeG_iContracted) :
            MinimalSpanningTree(G, storeOriginalGraph),
            fullBranchingTree(storeFullBranchingTree ? std::make_optional(FullBranchingTree(G)) : std::nullopt) {
        // we need to copy original graph in order to use "indexEdges" and not alter `originalGraph`
        Graph G_i(G, true, false);
        G_i.indexEdges(true);
        // NetworKit doesn't handle multiEdges, hence we don't check for their existence.
        G_i.removeSelfLoops();
        Bookkeeping bookKeeping(G_i, fullBranchingTree);
        count limitedSteps = storeG_iContracted ? 2 : std::numeric_limits<count>::max();
        /// `limitedSteps` condition used for KktMST
        while (G_i.numberOfNodes() > 1 && G_i.numberOfEdges() > 0 && limitedSteps-- > 0 ) {
            std::vector<WeightedUndirectedEdgeDecoupled::type> edgesForMst(G_i.upperNodeIdBound());

            /** Choose "the blue" edges **/
            G_i.forNodes([&](node u) {
                auto neighborRange = G_i.weightNeighborRange(u);
                auto lightestEdge = std::min_element(neighborRange.begin(), neighborRange.end(),
                                                     WeightedUndirectedEdgeDecoupled::less);
#ifdef MST_DEBUG
                assert((*lightestEdge).first != NetworKit::none && (*lightestEdge).first != u);
#endif
                edgesForMst[u] = {(*lightestEdge).first, (*lightestEdge).second};
            });

            // Add edges to MST
            for (index u = 0; u < G_i.upperNodeIdBound(); u++) {
                auto[v, weight] = edgesForMst[u];
                // origU and origV could be swapped - it's fine.
                // networkit bug - assertion fails (used to fail) sometimes! I think I've tackled it, but it's safer to
                // keep the assert.
                assert(G_i.edgeId(u, v) == G_i.edgeId(v, u));
                auto[origU, origV] = bookKeeping.getOriginalNodes(G_i.edgeId(u, v));
                if (!mst.hasEdge(origU, origV)) {
                    mst.addEdge(origU, origV, weight);
                    bookKeeping.addEdgeToF_i(u, v, weight);
                }
            }
            G_i = bookKeeping.prepareForNextIteration(G_i, edgesForMst);
        }
        if (storeG_iContracted) {
            // BoruvkaMST computed only partially. Probably used for KktMST.
            contractedG_iAfterLimitedSteps = {G_i};
            contractedG_iGetOriginalNodesFromEdgeId = {bookKeeping.E_iToE_0};
        }
    };
}
