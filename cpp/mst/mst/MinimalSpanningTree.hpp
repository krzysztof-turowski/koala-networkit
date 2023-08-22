#pragma once

#include "utils.hpp"
#include <networkit/graph/Graph.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <optional>

namespace MST {

    class MinimalSpanningTree {
    protected:
        std::optional<const Graph> originalGraph;
        Graph mst;
    public:
        [[nodiscard]]
        const Graph& getSpanningTree() const { return mst; }

        [[nodiscard]]
        const std::optional<const Graph>& getOriginalGraph() const { return originalGraph; }

        explicit MinimalSpanningTree(const Graph& G, bool storeOriginalGraph = false) :
                originalGraph(storeOriginalGraph ? std::make_optional(G) : std::nullopt),
                mst(G.upperNodeIdBound(), true, false) {
            if (G.isDirected()) {
                throw std::invalid_argument("MinimalSpanningTree requires undirected graph.");
            }
            if (!G.isWeighted()) {
                throw std::invalid_argument("MinimalSpanningTree requires weighted graph.");
            }
        }

        [[nodiscard]]
        std::vector<std::pair<NodePair, edgeweight>> getEdgesNotBelongingToMST() const {
            /// Edges might appear twice.
            if (!originalGraph.has_value()) {
                throw std::invalid_argument("Original graph is not provided.");
            }
            std::vector<std::pair<NodePair, edgeweight>> G_minus_M;
            std::vector<bool> visited(originalGraph->upperNodeIdBound(), false);
            originalGraph->forNodes([&](const node u) {
                mst.forNeighborsOf(u, [&](const node v, const edgeweight) {
                    visited[v] = true;
                });
                originalGraph->forNeighborsOf(u, [&](const node v, const edgeweight ew) {
                    if (!visited[v]) {
                        G_minus_M.push_back({{u, v}, ew});
                    }
                });
                mst.forNeighborsOf(u, [&](const node v, const edgeweight) {
                    visited[v] = false;
                });
            });
            return G_minus_M;
        }

        virtual ~MinimalSpanningTree() = default;
    };
}
