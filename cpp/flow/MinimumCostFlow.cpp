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
    return flow[e];
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

} /* namespace Koala */
