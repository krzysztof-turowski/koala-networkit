#include <flow/MinimumCostFlow.hpp>
#include <networkit/graph/Graph.hpp>

namespace Koala {

MinimumCostFlow::MinimumCostFlow(NetworKit::Graph &graph,
    edge_map<MCFEdgeParams>& ep, node_map<int>& np)
    : graph(graph), edge_params(ep), excess(np) {}

MinimumCostFlow::MinimumCostFlow(NetworKit::Graph &graph,
    edge_map<MCFEdgeParams>& ep)
    : graph(graph), edge_params(ep) {}

MinimumCostFlow::MinimumCostFlow(NetworKit::Graph &graph,
    edge_map<MCFEdgeParams>& ep, NetworKit::node s, NetworKit::node t, int fl)
    : graph(graph), edge_params(ep), excess({{s, fl}, {t, -fl}}) {}

int MinimumCostFlow::getFlow(const edge& e) {
    return flow[flowEdgeToId[e]];
}

int MinimumCostFlow::getMinCost() const {
    return min_cost;
}

bool MinimumCostFlow::isOk() const {
    return feasible;
}

void MinimumCostFlow::constructCirculation() {
    NetworKit::node n = graph.addNode();
    graph.forNodes(
        [&](NetworKit::node v) {
            if (excess[v] > 0) {
                graph.addEdge(n, v);
                edge_params[{n,v}].capacity = excess[v];
                edge_params[{v,n}].capacity = -excess[v];
            }
            if (excess[v] < 0) {
                graph.addEdge(v, n);
                edge_params[{n,v}].capacity = excess[v];
                edge_params[{v,n}].capacity = -excess[v];
            }
        });
}

void MinimumCostFlow::constructFlow() {
    // TODO
}

void MinimumCostFlow::makeUncapacitated() {
    if (!graph.isWeighted()) return;
    int nodes = graph.numberOfNodes();
    NetworKit::Graph newGraph(nodes, false, true, true);
    node_map<int> newB;
    edgeid_map<int> newCosts;

    NetworKit::node firstAdded = graph.upperNodeIdBound()+1;
    uncapacitatedNodesBounds = {firstAdded, firstAdded-1};
    graph.forEdges(
    [&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight weight, NetworKit::edgeid id) {
        auto k = newGraph.addNode();
        uncapacitatedNodesBounds.second = k;
        auto ukId = newGraph.upperEdgeIdBound();
        newGraph.addEdge(u,k);
        auto vkId = newGraph.upperEdgeIdBound();
        newGraph.addEdge(v,k);
        newB[u] = b[u];
        newB[v] = weight + b[v];
        newB[k] = -weight;

        newCosts[ukId] = 0;
        newCosts[vkId] = costs[id];
        uncapacitatedMapping[id] = vkId;
        flowEdgeToId[{u,v}] = vkId;
    });

    b = newB;
    graph = newGraph;
    costs = newCosts;
}

} /* namespace Koala */
