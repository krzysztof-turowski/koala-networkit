#pragma once

#include <map>
#include <utility>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace Koala {

struct MCFEdgeParams {
    int capacity;
    int cost;
};

using edge = std::pair<NetworKit::node, NetworKit::node>;

template<typename T>
using edge_map = std::map<Koala::edge, T>;

template<typename T> 
using edgeid_map = std::map<NetworKit::edgeid, T>;

template<typename T>
using node_map = std::map<NetworKit::node, T>;

class MinimumCostFlow : public NetworKit::Algorithm {
 public:
    
    MinimumCostFlow() {};
    // flow
    MinimumCostFlow(NetworKit::Graph&, edge_map<MCFEdgeParams>&,
        NetworKit::node, NetworKit::node, int);
    // flow
    MinimumCostFlow(NetworKit::Graph&, edge_map<MCFEdgeParams>&,
        node_map<int>&);
    // circulation
    MinimumCostFlow(NetworKit::Graph&, edge_map<MCFEdgeParams>&);

    void run() {
        hasRun = false;
        runImpl();
        hasRun = true;
    }
    bool isOk() const;

    int getFlow(const edge&);
    int getMinCost() const;
    
 protected:
    virtual void runImpl() = 0;
    void constructCirculation();
    void constructFlow();

    bool modifiedUncapacitated = false;
    edgeid_map<NetworKit::edgeid> uncapacitatedMapping;
    std::pair<NetworKit::node, NetworKit::node> uncapacitatedNodesBounds = {0,-1}; 

    // Makes graph uncapacitated, 
    // if the graph was uncapacitated before, then returns
    void makeUncapacitated();
    void makeConnected();

    NetworKit::Graph graph;
    edge_map<MCFEdgeParams> edge_params;
    node_map<int> b;
    node_map<int> excess;

    edge_map<NetworKit::edgeid> flowEdgeToId;
    edgeid_map<int> flow;
    edgeid_map<int> costs;
    int min_cost{0};
    bool feasible{true};
};

} /* namespace Koala */
