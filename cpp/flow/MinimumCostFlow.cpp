#include <flow/MinimumCostFlow.hpp>
#include <networkit/graph/Graph.hpp>


namespace Koala{

// MinimumCostFlow::MinimumCostFlow() : graph(NetworKit::Graph(0,false, true, false)){}
// MinimumCostFlow::MinimumCostFlow(NetworKit::Graph &graph) : graph(graph){}
MinimumCostFlow::MinimumCostFlow(NetworKit::Graph &graph, edge_map<MCFEdgeParams>& ep, node_map<int>& np) 
    : graph(graph), edge_params(ep), excess(np) {}

MinimumCostFlow::MinimumCostFlow(NetworKit::Graph &graph, edge_map<MCFEdgeParams>& ep) 
    : graph(graph), edge_params(ep) {}
    
MinimumCostFlow::MinimumCostFlow(NetworKit::Graph &graph, edge_map<MCFEdgeParams>& ep, NetworKit::node s, NetworKit::node t, int fl)
    : graph(graph), edge_params(ep), excess({{s,fl}, {t,-fl}}){} 


int MinimumCostFlow::getFlow(const edge& e){
    return flow[e];
}

int MinimumCostFlow::getMinCost() const{
    return min_cost;
} 

bool MinimumCostFlow::isOk() const{
    return feasible;
}

} /* namespace Koala */