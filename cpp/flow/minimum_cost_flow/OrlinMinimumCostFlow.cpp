#include <flow/minimum_cost_flow/OrlinMinimumCostFlow.hpp>

#include <algorithm>
#include <cassert>
#include <functional>
#include <limits>
#include <stack>
#include <unordered_map>
#include <utility>
#include <vector>

#include <flow/GoldbergTarjanPushRelabelMaximumFlow.hpp>
#include <structures/heap/FibonacciHeap.hpp>

using node = NetworKit::node;

namespace Koala {

int64_t OrlinMinimumCostFlow::getFlow(NetworKit::Edge const& edge) {
    if (maxflow.has_value()) {
        return maxflow->getFlow({edge.u, edge.v});
    }
    return 0;
}

std::unordered_map<NetworKit::Edge, int64_t> OrlinMinimumCostFlow::getMinCostFlow() const {
    return computed_flow;
}

void OrlinMinimumCostFlow::initialize() {
    uncapacitated_nodes_bounds.first = network.getGraph().upperNodeIdBound();
    network.makeUncapacitated();
    uncapacitated_nodes_bounds.second = network.getGraph().upperNodeIdBound();
    network.makeConnected();
    auto& graph = network.getGraph();
    nodes_number = graph.numberOfNodes();
    max_node_id = graph.upperNodeIdBound();
    original_graph = graph;
    excess.assign(max_node_id, 0);
    for (auto [key, value] : network.excess) {
        excess[key] = value;
    }
    potential.assign(max_node_id, 0);
    potential_computed = potential;
    distances.assign(max_node_id, {std::numeric_limits<int64_t>::max(), 0});
    visited.assign(max_node_id, false);
    edges.assign(2 * graph.numberOfEdges(), Edge());
    NetworKit::index edge_pair_index = 0;

    neighbors.assign(max_node_id, std::vector<NetworKit::index>());

    graph.forNodes([&](node u) {
        graph.forNeighborsOf(u, [&](node v) {
            node from = u;
            node to = v;
            int64_t cost = network.cost[{u, v}];
            int64_t capacity = network.capacity[{u, v}];
            neighbors[from].push_back(2*edge_pair_index);
            edges[2*edge_pair_index] = {
                from, to,
                cost, capacity, 0LL
            };

            neighbors[to].push_back(2*edge_pair_index + 1);
            edges[2*edge_pair_index + 1] = {
                to, from,
                -cost, 0LL, 0LL
            };
            ++edge_pair_index;
        });
    });
    original_edges = edges;
    make_reduced_costs_nonnegative();
}

bool OrlinMinimumCostFlow::is_added_uncapacitated(node v) const {
    auto [begin, end] = uncapacitated_nodes_bounds;
    return begin <= v && v < end;
}

void OrlinMinimumCostFlow::apply_potential() {
    for (Edge& edge : edges) {
        edge.cost += potential[edge.to] - potential[edge.from];
    }
    std::fill(potential.begin(), potential.end(), 0);
}

void OrlinMinimumCostFlow::contract_nodes(node u, node v) {
    if (u > v) std::swap(u, v);
    for (NetworKit::index edge_index : neighbors[v]) {
        Edge& edge = edges[edge_index];
        if (edge.to == u) {
            edge.capacity = 0;
            edge.flow = 0;
            edges[edge_index ^ 1].capacity = 0;
            edges[edge_index ^ 1].flow = 0;
        } else {
            edge.from = u;
            neighbors[u].push_back(edge_index);

            edges[edge_index ^ 1].to = u;
        }
    }

    excess[u] += excess[v];
    excess[v] = 0;
    --nodes_number;
    neighbors[v].clear();
}

void OrlinMinimumCostFlow::push_no_excess(NetworKit::index edge_index, int64_t amount) {
    edges[edge_index].flow += amount;
    edges[edge_index ^ 1].flow -= amount;
}

int64_t OrlinMinimumCostFlow::find_optimal_delta(int64_t delta) {
    for (const Edge& edge : edges) {
        if (edge.flow != 0) {
            return delta;
        }
    }
    int64_t new_delta = 1;
    for (int64_t node_excess : excess) {
        while (new_delta < node_excess) {
            new_delta <<= 1;
        }
    }
    return new_delta;
}

void OrlinMinimumCostFlow::contraction_phase(int64_t delta) {
    apply_potential();
    for (NetworKit::index edge_index = 0; edge_index < edges.size(); ++edge_index) {
        const Edge& edge = edges[edge_index];
        if (edge.from != edge.to
                && edge.flow >= 3 * delta * static_cast<int64_t>(nodes_number)) {
            contract_nodes(edge.from, edge.to);
            // { u, v, edge_cost }
            contracted_nodes.push({edge.from, edge.to, original_edges[edge_index].cost});
        }
    }
}

bool OrlinMinimumCostFlow::is_imbalanced() {
    for (const int64_t& node_excess : excess) {
        if (node_excess != 0) {
            return true;
        }
    }
    return false;
}

void OrlinMinimumCostFlow::uncontract_nodes_potential() {
    while (!contracted_nodes.empty()) {
        auto [u, v, cost] = contracted_nodes.top();
        contracted_nodes.pop();

        if (u < v) {
            potential_computed[v] = potential_computed[u] - cost;
        } else {
            potential_computed[u] = potential_computed[v] + cost;
        }
    }
}

void OrlinMinimumCostFlow::augmenting_phase(node s, node t, int64_t delta) {
    dijkstra(s, delta);
    node current_node = t;
    while (current_node != s) {
        NetworKit::index edge_index = distances[current_node].second;
        push_no_excess(edge_index, delta);
        current_node = edges[edge_index].from;
    }
    excess[s] -= delta;
    excess[t] += delta;

    for (node v = 0; v < max_node_id; v++) {
        if (distances[v].first != std::numeric_limits<int64_t>::max()) {
            potential[v] -= distances[v].first;
            potential_computed[v] -= distances[v].first;
        }
    }
}

void OrlinMinimumCostFlow::run_impl() {
    initialize();

    int64_t delta = std::numeric_limits<int64_t>::max();

    while (is_imbalanced()) {
        delta = std::min(delta, find_optimal_delta(delta));
        contraction_phase(delta);

        auto [uncapacitated_begin, uncapacitated_end] = uncapacitated_nodes_bounds;
        for (node u = uncapacitated_begin; u < uncapacitated_end; u++) {
            for (NetworKit::index edge_index : neighbors[u]) {
                NetworKit::index incoming_edge_index = edge_index ^ 1;
                const Edge& incoming_edge = edges[incoming_edge_index];
                node v = incoming_edge.from;
                while (excess[u] <= -ALPHA*delta && excess[v] >= ALPHA*delta
                       && incoming_edge.cost == 0
                       && incoming_edge.capacity >= incoming_edge.flow + delta) {
                    push_no_excess(incoming_edge_index, delta);
                    excess[v] -= delta;
                    excess[u] += delta;
                }
                while (excess[v] <= -ALPHA*delta && excess[u] >= ALPHA*delta
                       && edges[edge_index].cost == 0
                       && edges[edge_index].capacity >= edges[edge_index].flow + delta) {
                    push_no_excess(edge_index, delta);
                    excess[u] -= delta;
                    excess[v] += delta;
                }
            }
        }

        std::stack<node> excess_nodes, deficit_nodes;

        for (node v = 0; v < max_node_id; v++) {
            if (excess[v] >= ALPHA*delta) {
                excess_nodes.push(v);
            } else if (excess[v] <= -ALPHA*delta) {
                deficit_nodes.push(v);
            }
        }

        while (!excess_nodes.empty() && !deficit_nodes.empty()) {
            node s = excess_nodes.top();
            node t = deficit_nodes.top();

            augmenting_phase(s, t, delta);
            if (excess[s] < ALPHA*delta) {
                excess_nodes.pop();
            }
            if (excess[t] > -ALPHA*delta) {
                deficit_nodes.pop();
            }
        }

        delta /= 2;
    }

    uncontract_nodes_potential();
    compute_final_flows();
}

void OrlinMinimumCostFlow::dijkstra(node source, int64_t delta) {
    using HeapKey = std::pair<int64_t, node>;
    using Heap = FibonacciHeap<HeapKey, std::greater<HeapKey>>;

    Heap queue;
    std::vector<NetworKit::index> heap_handles(max_node_id, NetworKit::none);

    // Relaxes one residual arc out of u, returning its head if the distance improved.
    auto relax = [&](node u, NetworKit::index edge_index) -> node {
        const Edge& edge = edges[edge_index];
        if (edge.capacity < edge.flow + delta) {
            return NetworKit::none;
        }
        node v = edge.to;
        int64_t distance =
            distances[u].first + edge.cost - potential[u] + potential[v];
        if (distance >= distances[v].first) {
            return NetworKit::none;
        }
        distances[v] = {distance, edge_index};
        return v;
    };

    // Queues v at its current distance, or lowers its key in place if already queued.
    auto enqueue = [&](node v) {
        HeapKey key{distances[v].first, v};
        if (heap_handles[v] == NetworKit::none) {
            heap_handles[v] = *queue.push(key);
        } else {
            queue.update(Heap::iterator(heap_handles[v]), key);
        }
    };

    std::fill(
        distances.begin(), distances.end(),
        std::make_pair(std::numeric_limits<int64_t>::max(), NetworKit::index{0}));
    std::fill(visited.begin(), visited.end(), false);
    distances[source] = {0, 0};
    enqueue(source);

    while (!queue.empty()) {
        node u = queue.top().second;
        queue.pop();
        heap_handles[u] = NetworKit::none;

        if (visited[u]) continue;
        visited[u] = true;

        for (NetworKit::index edge_index : neighbors[u]) {
            node v = relax(u, edge_index);
            if (v == NetworKit::none) continue;

            if (!is_added_uncapacitated(v)) {
                enqueue(v);
                continue;
            }
            // Nodes added when splitting capacitated arcs are never settled on
            // their own: step straight through them to the next real node.
            for (NetworKit::index next_edge_index : neighbors[v]) {
                if (edges[next_edge_index].to == u) continue;
                node next_node = relax(v, next_edge_index);
                if (next_node != NetworKit::none) {
                    enqueue(next_node);
                }
            }
        }
    }
}

void OrlinMinimumCostFlow::make_reduced_costs_nonnegative() {
    constexpr int64_t INF = std::numeric_limits<int32_t>::max();
    std::vector<int64_t> initial_distances(max_node_id, INF);
    initial_distances[0] = 0;

    bool changed = true;
    for (NetworKit::count iteration = 0; iteration < nodes_number && changed; ++iteration) {
        changed = false;
        for (const Edge& edge : edges) {
            if (edge.from == edge.to || edge.capacity <= edge.flow) continue;
            if (initial_distances[edge.from] == INF) continue;
            int64_t length = edge.cost - potential[edge.from] + potential[edge.to];
            if (initial_distances[edge.from] + length < initial_distances[edge.to]) {
                initial_distances[edge.to] = initial_distances[edge.from] + length;
                changed = true;
            }
        }
    }

    for (node v = 0; v < max_node_id; ++v) {
        if (initial_distances[v] == INF) continue;
        potential[v] -= initial_distances[v];
        potential_computed[v] -= initial_distances[v];
    }
}

void OrlinMinimumCostFlow::compute_final_flows() {
    NetworKit::Graph maxflow_graph(max_node_id, true, true);
    original_graph.forEdges([&](node u, node v) {
        auto cost = network.cost[{u, v}];

        if (cost - potential_computed[u] + potential_computed[v] == 0) {
            NetworKit::edgeweight infinite_capacity = std::numeric_limits<std::int32_t>::max();
            maxflow_graph.addEdge(u, v, infinite_capacity);
        }
    });

    node s = maxflow_graph.addNode();
    node t = maxflow_graph.addNode();

    for (auto [key, value] : network.excess) {
        if (value > 0) {
            maxflow_graph.addEdge(s, key, value);
        } else if (value < 0) {
            maxflow_graph.addEdge(key, t, -value);
        }
    }
    maxflow.emplace(maxflow_graph, s, t);
    maxflow->run();
    min_cost = 0;
    computed_flow.clear();
    maxflow_graph.forEdges([&](node u, node v) {
        if (u == s || v == t) return;
        int64_t flow = maxflow->getFlow({u, v});
        computed_flow[{u, v}] = flow;
        min_cost += network.cost[{u, v}] * flow;
    });
}

} /* namespace Koala */
