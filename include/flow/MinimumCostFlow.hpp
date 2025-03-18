#pragma once

#include <map>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace Koala {

struct MinCostFlowEdgeParams{ 
    int capacity;
    int cost;
};

typedef std::map<std::pair<NetworKit::node, NetworKit::node>, MinCostFlowEdgeParams> edge_map;
typedef std::map<NetworKit::node, int> node_map;

struct MinCostFlowParameters{
    NetworKit::Graph graph;
    NetworKit::node s;
    NetworKit::node t;
    edge_map edge_params;
    node_map supply;
};

class MinimumCostFlow : public NetworKit::Algorithm {
public:
    MinimumCostFlow(MinCostFlowParameters params);
    MinimumCostFlow();
    
    void addEdge(const std::pair<NetworKit::node, NetworKit::node>&, MinCostFlowEdgeParams);
    MinCostFlowEdgeParams getEdgeParams(const std::pair<NetworKit::node, NetworKit::node>&) const;
    void setEdgeParams(const std::pair<NetworKit::node, NetworKit::node>&, MinCostFlowEdgeParams);

    int getSupply(const NetworKit::node&) const;
    void setSupply(const NetworKit::node&);

    virtual void run();
    int getFlow(const std::pair<NetworKit::node, NetworKit::node>&) const;

    int getMinCost() const;
};

}