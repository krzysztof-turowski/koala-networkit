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
 public:
    explicit MCFlowNetwork(NetworKit::Graph const& g, bool circulation = false);
    MCFlowNetwork(NetworKit::Graph const& g,
        std::unordered_map<NetworKit::Edge, int64_t> const& cost, bool circulation = false);
    MCFlowNetwork(
        NetworKit::Graph const& g, std::unordered_map<NetworKit::Edge, int64_t> const& cost,
        std::unordered_map<NetworKit::node, int64_t> const& node_excess, bool circulation = false);

    NetworKit::Graph& getGraph();
    NetworKit::node addNode(int64_t node_excess);
    void addEdge(NetworKit::node s, NetworKit::node t, int64_t cost = 0, int64_t capacity = 0);

    std::unordered_map<NetworKit::Edge, int64_t> cost;
    std::unordered_map<NetworKit::Edge, int64_t> capacity;
    std::unordered_map<NetworKit::node, int64_t> excess;

    void makeConnected();
    void makeUncapacitated();
    void makeCostsNonNegative();
};

} /* namespace Koala */
