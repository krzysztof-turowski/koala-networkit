#pragma once

#include <map>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace Koala {

struct MCFEdgeParams{ 
    int capacity;
    int cost;
};

template<typename T> 
using edge_map = std::map<std::pair<NetworKit::node, NetworKit::node>, T>;

typedef std::map<NetworKit::node, int> node_map;

class MinimumCostFlow : public NetworKit::Algorithm {
public:
    MinimumCostFlow();
    MinimumCostFlow(NetworKit::Graph&);
    MinimumCostFlow(NetworKit::Graph&, edge_map<MCFEdgeParams>&, node_map&);

    NetworKit::Graph& getGraph();
    MCFEdgeParams getEdgeParams(const std::pair<NetworKit::node, NetworKit::node>&);
    void setEdgeParams(const std::pair<NetworKit::node, NetworKit::node>&, MCFEdgeParams);

    int getSupply(const NetworKit::node&);
    void setSupply(const NetworKit::node&, int);

    virtual void run();
    int getFlow(const std::pair<NetworKit::node, NetworKit::node>&);

    int getMinCost() const;

protected:
    NetworKit::Graph graph;
    edge_map<MCFEdgeParams> edge_params;
    node_map supply;
    edge_map<int> flow;
    int min_cost{0};
};

}