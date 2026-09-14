#pragma once

#include <cstdint>
#include <unordered_map>
#include <vector>

#include <networkit/graph/EdgeIterators.hpp>
#include <networkit/graph/Graph.hpp>

namespace Koala {

class MCFlowNetwork {
    bool uncapacitated = false;
    NetworKit::Graph graph;
    std::unordered_map<NetworKit::Edge, NetworKit::Edge> uncapacitated_to_original_edge_mapping;
 public:
    explicit MCFlowNetwork(NetworKit::Graph const& g);
    MCFlowNetwork(NetworKit::Graph const& g,
        std::unordered_map<NetworKit::Edge, int64_t> const& cost);
    MCFlowNetwork(
        NetworKit::Graph const& g, std::unordered_map<NetworKit::Edge, int64_t> const& cost,
        std::unordered_map<NetworKit::node, int64_t> const& node_excess);

    NetworKit::Graph& getGraph();
    NetworKit::node addNode(int64_t node_excess);
    void addEdge(NetworKit::node s, NetworKit::node t, int64_t cost = 0, int64_t capacity = 0);

    std::unordered_map<NetworKit::Edge, int64_t> cost;
    std::unordered_map<NetworKit::Edge, int64_t> capacity;
    std::unordered_map<NetworKit::node, int64_t> excess;

    /**
     * Transforms the network to ensure that it is connected, which may involve adding edges or nodes as necessary.
     * This transformation is performed in-place, modifying the original network.
     * After this operation, the network will be connected, meaning there is a path between any two nodes in the network.
     * Note: This operation may change the structure of the network,
     */
    void makeConnected();
    /**
     * Transforms the network to ensure that it is uncapacitated, which may involve adding edges or nodes as necessary.
     * This transformation is performed in-place, modifying the original network.
     * After this operation, the network will be uncapacitated, meaning that all edges have infinite capacity.
     * Note: This operation may change the structure of the network,
     */
    void makeUncapacitated();
    /**
     * Transforms the network to ensure that all costs are non-negative, which may involve adding edges or nodes as necessary.
     * This transformation is performed in-place, modifying the original network.
     * After this operation, all costs in the network will be non-negative.
     * Note: This operation may change the structure of the network,
     */
    void makeCostsNonNegative();
    /** 
     * Returns the edge in the original network that corresponds 
     * to the given edge with associated cost in the transformed (uncapacitated) network .
     * @param edge The edge in the transformed (uncapacitated) network.
     * @return The corresponding edge with associated cost in the original network.
     */
    NetworKit::Edge getUncapacitatedToOriginalEdgeMapping(NetworKit::Edge const& edge) const;
};

} /* namespace Koala */
