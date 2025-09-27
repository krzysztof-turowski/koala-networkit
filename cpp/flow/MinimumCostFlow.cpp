#include <flow/MinimumCostFlow.hpp>
#include <networkit/graph/Graph.hpp>

namespace Koala {

using index = NetworKit::index;

MinimumCostFlow::MinimumCostFlow(
    NetworKit::Graph const& graph,
    edgeid_map<long long> const& costs,
    node_map<long long> const& b
) : graph(graph), costs(costs), b(b) {
}

MinimumCostFlow::MinimumCostFlow(
    NetworKit::Graph const& graph,
    edgeid_map<std::pair<long long,long long>> const& cap_bounds, 
    edgeid_map<long long> const& costs
) : graph(graph), costs(costs) { 
    for(auto [key, val] : cap_bounds) {
        lowerbound[key] = val.first;
        upperbound[key] = val.second;
    }
}

MinimumCostFlow::MinimumCostFlow(
    NetworKit::Graph const& graph,
    edgeid_map<long long> const& ep,
    NetworKit::node s, NetworKit::node t, 
    long long fl
) : graph(graph), costs(ep), b({{s, fl}, {t, -fl}}) {
}

long long MinimumCostFlow::getFlow(NetworKit::edgeid e) {
    return flow[e];
}

long long MinimumCostFlow::getFlow(edge const& e) {
    return flow[flowEdgeToId[e]];
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
                index id = graph.upperEdgeIdBound();
                graph.addEdge(n, v);
                upperbound[id] = excess[v];
                lowerbound[id] = excess[v];
            }
            if (b[v] < 0) {
                index id = graph.upperEdgeIdBound();
                graph.addEdge(v, n);
                upperbound[id] = excess[v];
                lowerbound[id] = excess[v];
            }
        });
}

void MinimumCostFlow::constructFlowFromCirculation() {
    // TODO    
    NetworKit::Graph G(graph.numberOfNodes(), true, true, true);

    graph.forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeid eid) {
        b[u] += lowerbound[eid];
        b[v] -= lowerbound[eid];
        NetworKit::index i = G.upperEdgeIdBound();
        G.addEdge(u, v, upperbound[eid] - lowerbound[eid]);
        circulationMapping[eid] = i;
    });

    graph = G;

}

void MinimumCostFlow::makeUncapacitated() {
    if (!graph.isWeighted()) return;
    int nodes = graph.numberOfNodes();
    NetworKit::Graph newGraph(nodes, false, true, true);
    node_map<long long> newB;
    edgeid_map<long long> newCosts;

    NetworKit::node firstAdded = graph.upperNodeIdBound();
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
        costs[bound] = costs[bound+1] = maxCost;
        graph.addEdge(u, sx, sumB);
        graph.addEdge(sx, u, sumB);
    });
}

} /* namespace Koala */
