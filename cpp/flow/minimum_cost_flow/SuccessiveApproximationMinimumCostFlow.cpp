#include <flow/minimum_cost_flow/SuccessiveApproximationMinimumCostFlow.hpp>

#include <algorithm>
#include <limits>
#include <memory>
#include <vector>

namespace Koala {

using node = NetworKit::node;

inline double SuccessiveApproximationMinimumCostFlow::reduced_cost(NetworKit::index edge_index) {
    const Edge& edge = edges[edge_index];
    return static_cast<double>(edge.cost) - potential[edge.from] + potential[edge.to];
}

inline int64_t SuccessiveApproximationMinimumCostFlow::residual_capacity(
        NetworKit::index edge_index) {
    return edges[edge_index].capacity - edges[edge_index].flow;
}

bool SuccessiveApproximationMinimumCostFlow::is_imbalanced() {
    for (int64_t node_excess : excess) {
        if (node_excess) return true;
    }
    return false;
}

void SuccessiveApproximationMinimumCostFlow::force_flow(
        NetworKit::index edge_index, int64_t amount) {
    Edge& edge = edges[edge_index];
    edge.flow += amount;
    edges[edge_index ^ 1].flow -= amount;
    excess[edge.from] -= amount;
    excess[edge.to] += amount;
}

void SuccessiveApproximationMinimumCostFlow::push(NetworKit::index edge_index) {
    node u = edges[edge_index].from;
    if (excess[u] > 0) {
        force_flow(edge_index, std::min(residual_capacity(edge_index), excess[u]));
    }
}

void SuccessiveApproximationMinimumCostFlow::relabel(NetworKit::node const& u) {
    double new_potential = std::numeric_limits<double>::infinity();

    for (NetworKit::index edge_index : neighbors[u]) {
        if (residual_capacity(edge_index) > 0) {
            const Edge& edge = edges[edge_index];
            new_potential = std::min(new_potential, potential[edge.to] + epsilon + edge.cost);
        }
    }

    potential[u] = new_potential;
}

void SuccessiveApproximationMinimumCostFlow::refine() {
    epsilon /= 2;
    for (NetworKit::index edge_index = 0; edge_index < edges.size(); ++edge_index) {
        if (reduced_cost(edge_index) < 0) {
            force_flow(edge_index, residual_capacity(edge_index));
        }
    }
    wave();
}

void SuccessiveApproximationMinimumCostFlow::wave() {
    std::unique_ptr<DischargeList> list = std::make_unique<ToposortList>(*this);

    node v = list->getNext();
    while (is_imbalanced()) {
        if (excess[v] > 0) {
            bool relabeled = discharge(v);
            if (relabeled) {
                list->moveToStart();
            }
        }
        v = list->getNext();
    }
}

bool SuccessiveApproximationMinimumCostFlow::discharge(NetworKit::node const& u) {
    int64_t& node_excess = excess[u];
    for (NetworKit::index edge_index : neighbors[u]) {
        if (node_excess && reduced_cost(edge_index) < 0
                && residual_capacity(edge_index) > 0) {
            push(edge_index);
        }
    }

    if (node_excess > 0) {
        relabel(u);
        return true;
    }

    return false;
}

void SuccessiveApproximationMinimumCostFlow::initialize() {
    auto& graph = network.getGraph();
    NetworKit::count max_node_id = graph.upperNodeIdBound();
    nodes_number = graph.numberOfNodes();
    potential.clear();
    excess.assign(max_node_id, 0);

    for (auto [key, value] : network.excess) {
        excess[key] = value;
    }

    potential.assign(max_node_id, 0);
    edges.reserve(2 * graph.numberOfEdges());

    neighbors.assign(max_node_id, std::vector<NetworKit::index>());
    int64_t max_cost = 1;

    graph.forNodes([&](node u) {
        graph.forNeighborsOf(u, [&](node v) {
            node from = u;
            node to = v;
            int64_t cost = network.cost[{u, v}];
            int64_t capacity = network.capacity[{u, v}];

            neighbors[from].push_back(edges.size());
            edges.push_back({
                from, to,
                cost, capacity, 0LL
            });

            max_cost = std::max(max_cost, std::abs(cost));

            neighbors[to].push_back(edges.size());
            edges.push_back({
                to, from,
                -cost, 0LL, 0LL
            });
        });
    });

    epsilon = static_cast<double>(max_cost);
}

void SuccessiveApproximationMinimumCostFlow::run_impl() {
    initialize();

    while (epsilon >= 1.0/nodes_number) {
        refine();
    }

    min_cost = 0;

    for (const Edge& edge : edges) {
        min_cost += edge.flow * edge.cost;
        computed_flow[{edge.from, edge.to}] = edge.flow;
    }
    min_cost /= 2;
}

SuccessiveApproximationMinimumCostFlow::ToposortList::ToposortList(
    SuccessiveApproximationMinimumCostFlow &algorithm) : algorithm(algorithm) {
    visited.assign(algorithm.nodes_number, 0);
    auto& graph = algorithm.network.getGraph();
    for (auto v : graph.nodeRange()) {
        if (!visited[v]) dfs(v);
    }
    next = nodes.begin();
}

void SuccessiveApproximationMinimumCostFlow::ToposortList::dfs(NetworKit::node u) {
    visited[u] = true;

    for (NetworKit::index edge_index : algorithm.neighbors[u]) {
        if (algorithm.reduced_cost(edge_index) < 0
                && algorithm.residual_capacity(edge_index) > 0) {
            node v = algorithm.edges[edge_index].to;
            if (!visited[v]) {
                dfs(v);
            }
        }
    }

    nodes.push_front(u);
}

NetworKit::node SuccessiveApproximationMinimumCostFlow::ToposortList::getNext() {
    if (next != nodes.end()) {
        current = next;
        next++;
    } else {
        next = nodes.begin();
        current = next++;
    }
    return *current;
}

void SuccessiveApproximationMinimumCostFlow::ToposortList::moveToStart() {
    if (current != nodes.end()) {
        nodes.splice(nodes.begin(), nodes, current);
    }
}

int64_t SuccessiveApproximationMinimumCostFlow::getFlow(NetworKit::Edge const& edge) {
    return computed_flow[edge];
}

} /* namespace Koala */
