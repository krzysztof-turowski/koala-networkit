#include <flow/MinimumCostFlow.hpp>
#include <networkit/graph/Graph.hpp>

namespace Koala {

using index = NetworKit::index;

MinimumCostFlow::MinimumCostFlow(
    NetworKit::Graph const& graph,
    edge_map<long long> const& costs,
    node_map<long long> const& b
) : graph(graph), costs(costs), b(b) {
}

MinimumCostFlow::MinimumCostFlow(
    NetworKit::Graph const& graph,
    edge_map<std::pair<long long,long long>> const& cap_bounds, 
    edge_map<long long> const& costs
) : graph(graph), costs(costs) { 
    for(auto [key, val] : cap_bounds) {
        lowerbound[key] = val.first;
        upperbound[key] = val.second;
    }
}

MinimumCostFlow::MinimumCostFlow(
    NetworKit::Graph const& graph,
    edge_map<long long> const& ep,
    NetworKit::node s, NetworKit::node t, 
    long long fl
) : graph(graph), costs(ep), b({{s, fl}, {t, -fl}}) {
}

long long MinimumCostFlow::getFlow(edge const& e) {
    edge ed = modifiedUncapacitated ? uncapacitatedMapping[e] : e;
    return flow[ed] + shouldAddLowerbounds ? lowerbound[ed] : 0;
}

long long MinimumCostFlow::getMinCost() const {
    return min_cost;
}

bool MinimumCostFlow::isOk() const {
    return feasible;
}

void MinimumCostFlow::constructCirculationFromFlow() {
    NetworKit::node n = graph.addNode();
    graph.forNodes(
        [&](NetworKit::node v) {
            if (b[v] > 0) {
                graph.addEdge(n, v);
                upperbound[{n, v}] = excess[v];
                lowerbound[{n, v}] = excess[v];
            }
            if (b[v] < 0) {
                graph.addEdge(v, n);
                upperbound[{v, n}] = excess[v];
                lowerbound[{v, n}] = excess[v];
            }
        });
}

void MinimumCostFlow::constructFlowFromCirculation() {
    NetworKit::Graph G(graph.numberOfNodes(), true, true, true);

    graph.forEdges([&](NetworKit::node u, NetworKit::node v) {
        b[u] -= lowerbound[{u, v}];
        b[v] += lowerbound[{u, v}];
        G.addEdge(u, v, upperbound[{u, v}] - lowerbound[{u, v}]);
    });

    graph = G;

}

void MinimumCostFlow::makeUncapacitated() {
    if (!graph.isWeighted()) return;
    int nodes = graph.numberOfNodes();
    NetworKit::Graph newGraph(nodes, false, true, true);
    node_map<long long> newB;
    edge_map<long long> newCosts;

    NetworKit::node firstAdded = graph.upperNodeIdBound();
    uncapacitatedNodesBounds = {firstAdded, firstAdded-1};
    graph.forEdges(
    [&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight weight) {
        auto k = newGraph.addNode();
        uncapacitatedNodesBounds.second = k;

        newGraph.addEdge(u,k);
        newGraph.addEdge(v,k);
        newB[u] = b[u];
        newB[v] = weight + b[v];
        newB[k] = -weight;

        newCosts[{u, k}] = 0;
        newCosts[{v, k}] = costs[{u, v}];
        uncapacitatedMapping[{u, v}] = {v, k};
        // flowEdgeToId[{u,v}] = {u, v};
    });

    b = newB;
    graph = newGraph;
    costs = newCosts;
}

void MinimumCostFlow::makeConnected() {
    long long maxCost{0}; 
    for (auto [edge, cost] : costs) {
        maxCost = std::max(maxCost, std::abs(cost));
    }
    
    maxCost *= graph.numberOfEdges() + 1;
    
    long long sumB{0};
    for (auto [nd, bval] : b) {
        if (bval > 0) sumB += bval;
    }

    NetworKit::node sx = graph.addNode();
    graph.forNodes([&](NetworKit::node u) {
        if (sx == u) return;
        auto bound = graph.upperEdgeIdBound();
        costs[{u, sx}] = costs[{sx, u}] = maxCost;
        graph.addEdge(u, sx, sumB);
        graph.addEdge(sx, u, sumB);
    });
}

} /* namespace Koala */
