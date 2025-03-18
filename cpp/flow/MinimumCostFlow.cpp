#include <flow/MinimumCostFlow.hpp>
#include <networkit/graph/Graph.hpp>

using edge = std::pair<NetworKit::node, NetworKit::node>;

namespace Koala{

MinimumCostFlow::MinimumCostFlow() : graph(NetworKit::Graph(0,false, true, false)){}
MinimumCostFlow::MinimumCostFlow(NetworKit::Graph &graph) : graph(graph){}
MinimumCostFlow::MinimumCostFlow(NetworKit::Graph &graph, edge_map<MCFEdgeParams>& ep, node_map& np) 
    : graph(graph), edge_params(ep), supply(np) {}

NetworKit::Graph& MinimumCostFlow::getGraph(){
    return graph;
}

MCFEdgeParams MinimumCostFlow::getEdgeParams(const edge &e){
    return edge_params[e];
}

void MinimumCostFlow::setEdgeParams(const edge& e, MCFEdgeParams p){
    edge_params[e] = p;
}

int MinimumCostFlow::getSupply(const NetworKit::node& n){
    return supply[n];
}

void MinimumCostFlow::setSupply(const NetworKit::node& n, int x){
    supply[n]=x;
}

int MinimumCostFlow::getFlow(const edge& e){
    return flow[e];
}

int MinimumCostFlow::getMinCost() const{
    return min_cost;
} 

} /* namespace Koala */