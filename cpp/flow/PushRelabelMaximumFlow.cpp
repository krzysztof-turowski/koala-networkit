#include <algorithm>
#include <limits>
#include <vector>

#include "flow/PushRelabelMaximumFlow.hpp"
#include "graph/GraphTools.hpp"

namespace Koala {

NetworKit::Edge reverse(const NetworKit::Edge &e) {
    return NetworKit::Edge(e.v, e.u);
}

void PushRelabelMaximumFlow::initialize() {
    capacity.clear(), flow.clear();
    distance.clear(), excess.clear();

    std::vector<NetworKit::WeightedEdge> original_edges;
    graph.forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight weight) {
        original_edges.emplace_back(u, v, weight);
    });
    for (const auto &e : original_edges) {
        capacity[e] = e.weight;
        if (!capacity.count(reverse(e))) {
            capacity[reverse(e)] = 0;
        }
    }
    for (const auto &e : original_edges) {
        if (!graph.hasEdge(e.v, e.u)) {
            graph.addEdge(e.v, e.u, 0);
        }
    }
    graph.forEdges([&](NetworKit::node u, NetworKit::node v) {
        flow[NetworKit::Edge(u, v)] = 0;
    });
    graph.forNodes([&](NetworKit::node v) {
        distance[v] = 0;
        excess[v] = 0;
    });
    distance[source] = graph.numberOfNodes();
    excess[source] = std::numeric_limits<int>::max();
}

void PushRelabelMaximumFlow::push(NetworKit::node u, NetworKit::node v) {
    const NetworKit::Edge e(u, v);
    const int delta = std::min(capacity[e] - flow[e], excess[u]);
    excess[u] -= delta, excess[v] += delta;
    flow[e] += delta, flow[reverse(e)] -= delta;
}

void PushRelabelMaximumFlow::relabel(NetworKit::node v) {
    int minimum = std::numeric_limits<int>::max();
    graph.forNeighborsOf(v, [&](NetworKit::node u) {
        const NetworKit::Edge e(v, u);
        if (capacity[e] - flow[e] > 0) {
            minimum = std::min(minimum, distance[u]);
        }
    });
    if (minimum != std::numeric_limits<int>::max()) {
        const int old_distance = distance[v];
        distance[v] = minimum + 1;
        if (distance[v] != old_distance) {
            on_relabel(v, old_distance);
        }
    }
}

void PushRelabelMaximumFlow::run() {
    GraphTools::ensureDirectedGraph(graph);
    initialize();
    graph.forNeighborsOf(source, [&](NetworKit::node v) {
        if (capacity[NetworKit::Edge(source, v)] > 0) {
            push(source, v);
        }
    });
    for (auto v = get_active_vertex(); v != NetworKit::none; v = get_active_vertex()) {
        const auto u = get_admissible_residual_edge(v);
        if (u == NetworKit::none) {
            relabel(v);
        } else {
            push(v, u);
        }
    }
    flow_size = excess[target];
    hasRun = true;
}

}  // namespace Koala
