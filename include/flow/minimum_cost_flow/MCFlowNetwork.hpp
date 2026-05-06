#pragma once

#include <networkit/graph/Graph.hpp>
#include <networkit/graph/EdgeIterators.hpp>
#include <unordered_map>
#include <vector>

namespace Koala {

class MCFlowNetwork {
    bool uncapacitated = false;
    NetworKit::Graph graph;
 public:
    MCFlowNetwork(NetworKit::Graph const& g, bool circulation = false);
    MCFlowNetwork(NetworKit::Graph const& g,
        std::unordered_map<NetworKit::Edge, int64_t> const& cost, bool circulation = false);
    MCFlowNetwork(NetworKit::Graph const& g, std::unordered_map<NetworKit::Edge, int64_t> const& cost, 
        std::unordered_map<NetworKit::node, int64_t> const& ex, bool circulation = false);

    NetworKit::Graph& getGraph();
    NetworKit::node addNode(std::int64_t ex);
    void addEdge(NetworKit::node s, NetworKit::node t, std::int64_t cost = 0, std::int64_t capacity = 0);

    std::unordered_map<NetworKit::Edge, std::int64_t> cost;
    std::unordered_map<NetworKit::Edge, std::int64_t> capacity;
    std::unordered_map<NetworKit::node, std::int64_t> excess;

    void makeConnected();
    void makeUncapacitated();
    void makeCostsNonNegative();
};

} /* namespace Koala */