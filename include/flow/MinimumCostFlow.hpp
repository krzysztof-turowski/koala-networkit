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
        run_impl();
        hasRun = true;
    }
    bool isOk() const;

    long long getFlow(edge const&);
    long long getMinCost() const;
    
 protected:
    virtual void run_impl() = 0;
    virtual void init_circulation() = 0;
    virtual void init_flow() = 0;

    void construct_circulation_from_flow();
    void construct_flow_from_circulation();

    bool modified_uncapacitated = false;
    edge_map<std::pair<NetworKit::node, NetworKit::node>> uncapacitated_mapping;
    std::pair<NetworKit::node, NetworKit::node> uncapacitated_nodes_bounds = {0,-1}; 
    bool is_added_uncapacitated(NetworKit::node);
    bool should_add_lowerbounds = false;
    // bool fromCirculation = false;
    // edge_map<NetworKit::edgeid> circulationMapping;

    // Makes graph uncapacitated, 
    // if the graph was uncapacitated before, then returns
    void make_uncapacitated();
    void make_connected();

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
