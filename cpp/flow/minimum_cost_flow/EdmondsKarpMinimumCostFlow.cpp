#include <flow/minimum_cost_flow/EdmondsKarpMinimumCostFlow.hpp>

#include <cassert>
#include <functional>
#include <limits>
#include <stack>
#include <unordered_map>
#include <utility>
#include <vector>

#include <structures/heap/FibonacciHeap.hpp>

using node = NetworKit::node;

namespace Koala {

void EdmondsKarpMinimumCostFlow::initialize() {
    auto& graph = network.getGraph();
    network.makeConnected();
    max_node_id = graph.upperNodeIdBound();

    excess.assign(max_node_id, 0);
    for (auto [key, value] : network.excess) {
        excess[key] = value;
    }
    potential.assign(max_node_id, 0);
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
}

std::vector<std::pair<int64_t, NetworKit::index>> EdmondsKarpMinimumCostFlow::dijkstra(
        node source, int64_t delta) {
    using HeapKey = std::pair<int64_t, node>;
    using Heap = FibonacciHeap<HeapKey, std::greater<HeapKey>>;

    Heap queue;
    std::vector<NetworKit::index> heap_handles(max_node_id, NetworKit::none);
    std::vector<std::pair<int64_t, NetworKit::index>> distances(
        max_node_id, {std::numeric_limits<int64_t>::max(), 0});
    std::vector<bool> visited(max_node_id, false);

    // Queues v at its current distance, or lowers its key in place if already queued.
    auto enqueue = [&](node v) {
        HeapKey key{distances[v].first, v};
        if (heap_handles[v] == NetworKit::none) {
            heap_handles[v] = *queue.push(key);
        } else {
            queue.update(Heap::iterator(heap_handles[v]), key);
        }
    };

    distances[source] = {0, 0};
    enqueue(source);

    while (!queue.empty()) {
        auto [distance, u] = queue.top();
        queue.pop();
        heap_handles[u] = NetworKit::none;

        if (visited[u]) continue;
        visited[u] = true;

        for (NetworKit::index edge_index : neighbors[u]) {
            const Edge& edge = edges[edge_index];
            if (edge.capacity >= edge.flow + delta) {
                assert(edge.from == u);
                int64_t new_distance = distance + edge.cost - potential[u] + potential[edge.to];
                if (new_distance < distances[edge.to].first) {
                    distances[edge.to] = {new_distance, edge_index};
                    enqueue(edge.to);
                }
            }
        }
    }

    return distances;
}

void EdmondsKarpMinimumCostFlow::send(NetworKit::index edge_index, int64_t amount) {
    edges[edge_index].flow += amount;
    edges[edge_index ^ 1].flow -= amount;
}

void EdmondsKarpMinimumCostFlow::augmenting_phase(node s, node t, int64_t delta) {
    auto distances = dijkstra(s, delta);
    node current_node = t;
    while (s != current_node) {
        auto [_, edge_index] = distances[current_node];
        send(edge_index, delta);
        current_node = edges[edge_index].from;
    }
    excess[t] += delta;
    excess[s] -= delta;

    for (node v = 0; v < max_node_id; v++) {
        if (distances[v].first != std::numeric_limits<int64_t>::max()) {
            potential[v] -= distances[v].first;
        }
    }
}

void EdmondsKarpMinimumCostFlow::delta_scaling_phase(int64_t delta) {
    NetworKit::Graph const& graph = network.getGraph();
    for (NetworKit::index edge_index = 0; edge_index < edges.size(); ++edge_index) {
        Edge& edge = edges[edge_index];
        if (edge.capacity - edge.flow >= delta
            && edge.cost - potential[edge.from] + potential[edge.to] <= 0) {
                send(edge_index, delta);
                excess[edge.from] -= delta;
                excess[edge.to] += delta;
        }
    }

    std::stack<node> excess_nodes, deficit_nodes;
    graph.forNodes([&](node v){
        int64_t node_excess = excess[v];
        if (node_excess >= delta) excess_nodes.push(v);
        else if (node_excess <= -delta) deficit_nodes.push(v);
    });

    while (!excess_nodes.empty() && !deficit_nodes.empty()) {
        node s = excess_nodes.top();
        node t = deficit_nodes.top();

        augmenting_phase(s, t, delta);
        // update excess_nodes, deficit_nodes
        if (excess[s] < delta) excess_nodes.pop();
        if (excess[t] > -delta) deficit_nodes.pop();
    }
}

void EdmondsKarpMinimumCostFlow::run_impl() {
    initialize();

    int64_t delta{1};
    for (int64_t node_excess : excess) {
        while (delta < node_excess) {
            delta <<= 1;
        }
    }

    while (delta >= 1) {
        delta_scaling_phase(delta);
        delta /= 2;
    }

    min_cost = 0;
    for (const Edge& edge : edges) {
        computed_flow[{edge.from, edge.to}] = edge.flow;
        min_cost += edge.flow * edge.cost;
    }
    min_cost /= 2;
}

std::unordered_map<NetworKit::Edge, int64_t> EdmondsKarpMinimumCostFlow::getMinCostFlow() const {
    return computed_flow;
}

int64_t EdmondsKarpMinimumCostFlow::getFlow(NetworKit::Edge const& edge) {
    return computed_flow[edge];
}

} /* namespace Koala */
