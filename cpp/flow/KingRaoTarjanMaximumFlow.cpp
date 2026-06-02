#include <utility>
#include <vector>

#include "flow/KingRaoTarjanMaximumFlow.hpp"
#include "graph/GraphTools.hpp"

namespace Koala {

KingRaoTarjanMaximumFlow::KingRaoTarjanMaximumFlow(
    NetworKit::Graph &graph, NetworKit::node source, NetworKit::node target,
    KRTEdgeDesignator::Parameters edge_designator_parameters)
    : PushRelabelMaximumFlow(graph, source, target),
      edge_designator_parameters(std::move(edge_designator_parameters)) { }

void KingRaoTarjanMaximumFlow::run() {
    GraphTools::ensureDirectedGraph(graph);
    NetworKit::Graph residual_graph(graph);
    std::vector<NetworKit::Edge> edges;
    residual_graph.forEdges([&](NetworKit::node u, NetworKit::node v) {
        edges.emplace_back(u, v);
    });
    for (const auto &edge : edges) {
        if (!residual_graph.hasEdge(edge.v, edge.u)) {
            residual_graph.addEdge(edge.v, edge.u);
        }
    }
    edge_designator.initialize(residual_graph, edge_designator_parameters);
    PushRelabelMaximumFlow::run();
}

void KingRaoTarjanMaximumFlow::on_relabel(NetworKit::node v, int old_distance) {
    edge_designator.response_adversary(v, old_distance);
}

NetworKit::node KingRaoTarjanMaximumFlow::get_active_vertex() {
    for (NetworKit::node v = 0; v < graph.upperNodeIdBound(); ++v) {
        if (v != source && v != target && graph.hasNode(v) && excess[v] > 0) {
            return v;
        }
    }
    return NetworKit::none;
}

NetworKit::node KingRaoTarjanMaximumFlow::get_admissible_residual_edge(NetworKit::node v) {
    for (auto u = edge_designator.current_edge(v, distance[v]); u != NetworKit::none;
            u = edge_designator.current_edge(v, distance[v])) {
        const NetworKit::Edge e(v, u);
        if (capacity[e] - flow[e] > 0 && distance[v] == distance[u] + 1) {
            return u;
        }
        edge_designator.response_adversary(v, distance[v], u, distance[v] - 1);
    }
    return NetworKit::none;
}

}  // namespace Koala
