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
    
    // MinimumCostFlow() {};
    // flow
    MinimumCostFlow(NetworKit::Graph const&, edge_map<long long> const&,
        NetworKit::node, NetworKit::node, long long);
    // flow
    MinimumCostFlow(NetworKit::Graph const&, edge_map<long long> const&,
        node_map<long long> const&);
    // circulation
    MinimumCostFlow(NetworKit::Graph const&, edge_map<std::pair<long long,long long>> const&,
        edge_map<long long> const&);

    void run() {
        hasRun = false;
        runImpl();
        hasRun = true;
    }
    bool isOk() const;

    long long getFlow(edge const&);
    long long getMinCost() const;
    
 protected:
    virtual void runImpl() = 0;
    virtual void initCirculation() = 0;
    virtual void initFlow() = 0;
    void initFlowWrapper() { initFlow();}
    void initCirculationWrapper() { initFlow();}

    void constructCirculationFromFlow();
    void constructFlowFromCirculation();

    bool modifiedUncapacitated = false;
    edge_map<std::pair<NetworKit::node, NetworKit::node>> uncapacitatedMapping;
    std::pair<NetworKit::node, NetworKit::node> uncapacitatedNodesBounds = {0,-1}; 

    bool shouldAddLowerbounds = false;
    // bool fromCirculation = false;
    // edge_map<NetworKit::edgeid> circulationMapping;

    // Makes graph uncapacitated, 
    // if the graph was uncapacitated before, then returns
    void makeUncapacitated();
    void makeConnected();

    NetworKit::Graph graph;
    edge_map<MCFEdgeParams> edge_params;
    // flow
    node_map<long long> b;
    node_map<long long> excess;
    node_map<long long> potential; 

    edge_map<NetworKit::edgeid> flowEdgeToId;
    edge_map<long long> flow;
    edge_map<long long> costs;
    edge_map<long long> upperbound;
    edge_map<long long> lowerbound;

    long long min_cost{0};
    bool feasible{true};
};

} /* namespace Koala */
