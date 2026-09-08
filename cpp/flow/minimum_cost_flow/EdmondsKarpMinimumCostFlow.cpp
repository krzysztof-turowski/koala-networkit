#include <flow/minimum_cost_flow/EdmondsKarpMinimumCostFlow.hpp>

#include <climits>
#include <set>
#include <stack>
#include <vector>

using node = NetworKit::node;
using edgeid = NetworKit::edgeid;
using edgeweight = NetworKit::edgeweight;

namespace Koala {

void EdmondsKarpMinimumCostFlow::initialize() {
    auto& graph = network.getGraph();
    network.makeConnected();
    n = graph.upperNodeIdBound();

    excess.assign(n, 0);
    for (auto [key, value] : network.excess) {
        excess[key] = value;
    }
    potential.assign(n, 0);
    edges.assign(2 * graph.numberOfEdges(), Edge());
    NetworKit::index ptr = 0;

    neighbors.assign(n, std::vector<NetworKit::index>());
    graph.forNodes([&](node u) {
        graph.forNeighborsOf(u, [&](node v) {
            NetworKit::node from = u;
            NetworKit::node to = v;
            int64_t cost = network.cost[{u, v}];
            int64_t capacity = network.capacity[{u, v}];
            neighbors[from].push_back(2*ptr);
            edges[2*ptr] = {
                from, to,
                cost, capacity, 0LL
            };

            neighbors[to].push_back(2*ptr+1);
            edges[2*ptr + 1] = {
                to, from,
                -cost, 0LL, 0LL
            };
            ++ptr;
        });
    });
}

std::vector<std::pair<int64_t, NetworKit::index>> EdmondsKarpMinimumCostFlow::dijkstra(
        node source, int64_t delta) {
    std::set<std::pair<int64_t, node>> pq;
    std::vector<std::pair<int64_t, NetworKit::index>> dist(
        n, {std::numeric_limits<int64_t>::max(), 0});
    dist[source] = {0, 0};
    pq.insert({0, source});
    std::vector<bool> visited(n, false);
    while (!pq.empty()) {
        auto [d, u] = *pq.begin();
        pq.erase(pq.begin());

        if (visited[u]) continue;
        visited[u] = true;

        for (NetworKit::index edge_idx : neighbors[u]) {
            const Edge& edge = edges[edge_idx];
            if (edge.capacity >= edge.flow + delta) {
                assert(edge.from == u);
                int64_t new_dist = d + edge.cost - potential[u] + potential[edge.to];
                if (new_dist < dist[edge.to].first) {
                    pq.erase({dist[edge.to].first, edge.to});
                    dist[edge.to] = {new_dist, edge_idx};
                    pq.insert({new_dist, edge.to});
                }
            }
        }
    }

    return dist;
}

void EdmondsKarpMinimumCostFlow::send(NetworKit::index edge_idx, int64_t value) {
    edges[edge_idx].flow += value;
    edges[edge_idx^1].flow -= value;
}

void EdmondsKarpMinimumCostFlow::augmenting_phase(node s, node t, int64_t delta) {
    auto dist = dijkstra(s, delta);
    node ptr = t;
    while (s != ptr) {
        auto [_, edge_idx] = dist[ptr];
        send(edge_idx, delta);
        ptr = edges[edge_idx].from;
    }
    excess[t] += delta;
    excess[s] -= delta;

    for (node i = 0; i < n; i++) {
        if (dist[i].first != LLONG_MAX) {
            potential[i] -= dist[i].first;
        }
    }
}

void EdmondsKarpMinimumCostFlow::delta_scaling_phase(int64_t delta) {
    NetworKit::Graph const& graph = network.getGraph();
    for (NetworKit::index i = 0; i < edges.size(); ++i) {
        Edge& edge = edges[i];
        if (edge.capacity - edge.flow >= delta
            && edge.cost - potential[edge.from] + potential[edge.to] <= 0) {
                send(i, delta);
                excess[edge.from] -= delta;
                excess[edge.to] += delta;
        }
    }

    std::stack<node> S, T;
    graph.forNodes([&](node v){
        int64_t ex = excess[v];
        if (ex >= delta) S.push(v);
        else if (ex <= -delta) T.push(v);
    });

    while (!S.empty() && !T.empty()) {
        node k = S.top();
        node l = T.top();

        augmenting_phase(k, l, delta);
        // update S, T
        if (excess[k] < delta) S.pop();
        if (excess[l] > -delta) T.pop();
    }
}

void EdmondsKarpMinimumCostFlow::run_impl() {
    initialize();

    int64_t delta{1};
    for (auto e : excess) {
        while (delta < e) {
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

int64_t EdmondsKarpMinimumCostFlow::getFlow(NetworKit::Edge const& edge) {
    return computed_flow[edge];
}

} /* namespace Koala */
