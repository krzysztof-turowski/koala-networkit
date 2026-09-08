#include <flow/minimum_cost_flow/MCFlowNetwork.hpp>
#include <networkit/graph/GraphTools.hpp>

#include <algorithm>
#include <limits>
#include <unordered_map>

using node = NetworKit::node;
using Edge = NetworKit::Edge;
using Graph = NetworKit:: Graph;
using int64 = std::int64_t;

namespace Koala {

MCFlowNetwork::MCFlowNetwork(Graph const& g, bool circulation = false) : graph(g) {
    if (graph.isWeighted()) {
        graph.forEdges([&](node u, node v, NetworKit::edgeweight weight, NetworKit::edgeid _) {
            capacity[{u, v}] += static_cast<int64>(weight < 0 ? weight - 0.5 : weight + 0.5);
        });
    } else {
        graph.forEdges([&](node u, node v) {
            capacity[{u, v}] = std::numeric_limits<int64>::max();
        });
    }
}

MCFlowNetwork::MCFlowNetwork(
    Graph const& g, std::unordered_map<Edge, int64> const& cost, bool circulation = false)
    : MCFlowNetwork(g, circulation) {
    this->cost = cost;
}

MCFlowNetwork::MCFlowNetwork(
    Graph const& g, std::unordered_map<Edge, int64> const& cost,
    std::unordered_map<node, int64> const& ex, bool circulation = false)
    : MCFlowNetwork(g, cost, circulation) {
    excess = ex;
}

Graph& MCFlowNetwork::getGraph() {
    return graph;
}

node MCFlowNetwork::addNode(int64 ex = 0) {
    node newNode = graph.addNode();
    excess[newNode] = ex;
    return newNode;
}

void MCFlowNetwork::addEdge(node s, node t, int64 cost = 0, int64 capacity = 0) {
    graph.addEdge(s, t, capacity);
    if (graph.isWeighted())
        this->capacity[{s, t}] += capacity;
    this->cost[{s, t}] = cost;
}

void MCFlowNetwork::makeConnected() {
    int64 max_cost{0};
    for (auto [edge, cost] : cost) {
        max_cost = std::max(max_cost, (int64)std::abs(cost));
    }

    max_cost *=  graph.numberOfEdges() + 1;

    NetworKit::node sx = graph.addNode();
    NetworKit::node sx2 = graph.addNode();
    graph.addEdge(sx, sx2, std::numeric_limits<NetworKit::edgeweight>::max());
    capacity[{sx, sx2}] = std::numeric_limits<int64>::max();
    cost[{sx, sx2}] = max_cost;
    graph.forNodes([&](NetworKit::node u) {
        if (sx == u || sx2 == u) return;
        auto bound = graph.upperEdgeIdBound();
        this->capacity[{u, sx}] = this->capacity[{sx2, u}] = std::numeric_limits<int64>::max();
        graph.addEdge(u, sx, std::numeric_limits<NetworKit::edgeweight>::infinity());
        graph.addEdge(sx2, u, std::numeric_limits<NetworKit::edgeweight>::infinity());
    });
}

void MCFlowNetwork::makeUncapacitated() {
    NetworKit::Graph g = NetworKit::GraphTools::copyNodes(graph);
    

    graph.forEdges([&](node u, node v, NetworKit::edgeweight weight) {
        int64 cap = capacity[{u, v}];
        if (cap <= 0 || cap == std::numeric_limits<int64>::max())
            return;
        node w = g.addNode();
        excess[v] += cap;
        excess[w] -= cap;
        capacity.erase({u, v});
        capacity[{u, w}] = capacity[{v, w}] = std::numeric_limits<int64>::max();
        cost[{u, w}] = cost[{u, v}];
        cost.erase({u, v});
        g.addEdge(u, w, std::numeric_limits<NetworKit::edgeweight>::infinity());
        g.addEdge(v, w, std::numeric_limits<NetworKit::edgeweight>::infinity());
    });

    graph = g;
    uncapacitated = true;
}

void MCFlowNetwork::makeCostsNonNegative() {
    NetworKit::Graph g(graph.upperNodeIdBound(), true, true);

    graph.forEdges([&](node u, node v) {
        int64 cap = capacity[{u, v}];
        if (cost[{u, v}] < 0) {
            excess[v] += cap;
            excess[u] -= cap;
            g.addEdge(v, u, static_cast<NetworKit::edgeweight>(cap));
            capacity[{v, u}] = capacity[{u, v}];
            capacity[{u, v}] = 0;
            cost[{v, u}] = -cost[{u, v}];
            cost[{u, v}] = 0;
        } else {
            g.addEdge(u, v, static_cast<NetworKit::edgeweight>(cap));
        }
    });

    graph = g;
}

} /* namespace Koala */
