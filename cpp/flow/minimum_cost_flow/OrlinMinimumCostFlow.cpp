#include <flow/minimum_cost_flow/OrlinMinimumCostFlow.hpp>

#include <algorithm>
#include <functional>
#include <limits>
#include <queue>
#include <stack>
#include <utility>
#include <vector>

#include <flow/GoldbergTarjanPushRelabelMaximumFlow.hpp>

using edgeid = NetworKit::edgeid;
using node = NetworKit::node;
using Edge = NetworKit::Edge;
using int64 = std::int64_t;

namespace Koala {

int64 OrlinMinimumCostFlow::getFlow(const NetworKit::Edge& edge) {
    if (maxflow.has_value()) {
        return maxflow->getFlow({edge.u, edge.v});
    }
    return 0;
}

void OrlinMinimumCostFlow::initialize() {
    uncapacitated_nodes_bounds.first = network.getGraph().upperNodeIdBound();
    network.makeUncapacitated();
    uncapacitated_nodes_bounds.second = network.getGraph().upperNodeIdBound();
    network.makeConnected();
    auto& graph = network.getGraph();
    nodes_number = graph.numberOfNodes();
    max_nodeid = graph.upperNodeIdBound();
    original_graph = graph;
    excess.assign(max_nodeid, 0);
    for (auto [key, value] : network.excess) {
        excess[key] = value;
    }
    potential.assign(max_nodeid, 0);
    potential_computed = potential;
    dist.assign(max_nodeid, {std::numeric_limits<int64_t>::max(), 0});
    visited.assign(max_nodeid, false);
    edges.assign(2 * graph.numberOfEdges(), Edge());
    NetworKit::index ptr = 0;

    neighbors.assign(max_nodeid, std::vector<NetworKit::index>());

    graph.forNodes([&](node u) {
        graph.forNeighborsOf(u, [&](node v) {
            node from = u;
            node to = v;
            int64 cost = network.cost[{u, v}];
            int64 capacity = network.capacity[{u, v}];
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
    for (NetworKit::index edge_idx : neighbors[v]) {
        Edge& edge = edges[edge_idx];
        if (edge.to == u) {
            edge.capacity = 0;
            edge.flow = 0;
            edges[edge_idx ^ 1].capacity = 0;
            edges[edge_idx ^ 1].flow = 0;
        } else {
            edge.from = u;
            neighbors[u].push_back(edge_idx);

            edges[edge_idx ^ 1].to = u;
        }
    }

    excess[u] += excess[v];
    excess[v] = 0;
    --nodes_number;
    neighbors[v].clear();
}

void OrlinMinimumCostFlow::push_no_excess(NetworKit::index edge_idx, int64_t amount) {
    edges[edge_idx].flow += amount;
    edges[edge_idx ^ 1].flow -= amount;
}

int64_t OrlinMinimumCostFlow::find_optimal_delta(int64_t delta) {
    for (const Edge& edge : edges) {
        if (edge.flow != 0) {
            return delta;
        }
    }
    int64_t newDelta = 1;
    for (int64_t excess : excess) {
        while (newDelta < excess) {
            newDelta <<= 1;
        }
    }
    return newDelta;
}

void OrlinMinimumCostFlow::contraction_phase(int64_t delta) {
    apply_potential();
    for (NetworKit::index i = 0; i < edges.size(); ++i) {
        const Edge& edge = edges[i];
        if (edge.from != edge.to
                && edge.flow >= 3 * delta * static_cast<int64_t>(nodes_number)) {
            contract_nodes(edge.from, edge.to);
            // { u, v, edge_cost }
            contracted_nodes.push({edge.from, edge.to, original_edges[i].cost});
        }
    }
}

bool OrlinMinimumCostFlow::is_imbalanced() {
    for (const int64_t& supply : excess) {
        if (supply != 0) {
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
    node ptr = t;
    while (ptr != s) {
        NetworKit::index edge_idx = dist[ptr].second;
        push_no_excess(edge_idx, delta);
        ptr = edges[edge_idx].from;
    }
    excess[s] -= delta;
    excess[t] += delta;

    for (node i = 0; i < max_nodeid; i++) {
        if (dist[i].first != std::numeric_limits<int64_t>::max()) {
            potential[i] -= dist[i].first;
            potential_computed[i] -= dist[i].first;
        }
    }
}

void OrlinMinimumCostFlow::run_impl() {
    initialize();

    int64_t delta = std::numeric_limits<int64_t>::max();

    while (is_imbalanced()) {
        delta = std::min(delta, find_optimal_delta(delta));
        contraction_phase(delta);

        auto [uncap_begin, uncap_end] = uncapacitated_nodes_bounds;
        for (node u = uncap_begin; u < uncap_end; u++) {
            for (NetworKit::index edge_idx : neighbors[u]) {
                NetworKit::index in_arc = edge_idx ^ 1;
                const Edge& e = edges[in_arc];
                node v = e.from;
                while (excess[u] <= -ALPHA*delta && excess[v] >= ALPHA*delta
                       && e.cost == 0
                       && e.capacity >= e.flow + delta) {
                    push_no_excess(in_arc, delta);
                    excess[v] -= delta;
                    excess[u] += delta;
                }
                while (excess[v] <= -ALPHA*delta && excess[u] >= ALPHA*delta
                       && edges[edge_idx].cost == 0
                       && edges[edge_idx].capacity >= edges[edge_idx].flow + delta) {
                    push_no_excess(edge_idx, delta);
                    excess[u] -= delta;
                    excess[v] += delta;
                }
            }
        }

        std::stack<node> S, T;

        for (node i = 0; i < max_nodeid; i++) {
            if (excess[i] >= ALPHA*delta) {
                S.push(i);
            } else if (excess[i] <= -ALPHA*delta) {
                T.push(i);
            }
        }

        while (!S.empty() && !T.empty()) {
            node s = S.top();
            node t = T.top();

            augmenting_phase(s, t, delta);
            if (excess[s] < ALPHA*delta) {
                S.pop();
            }
            if (excess[t] > -ALPHA*delta) {
                T.pop();
            }
        }

        delta /= 2;
    }

    uncontract_nodes_potential();
    compute_final_flows();
}

void OrlinMinimumCostFlow::dijkstra(node source, int64_t delta) {
    std::priority_queue<std::pair<int64_t, node>,
                        std::vector<std::pair<int64_t, node>>,
                        std::greater<>> pq;

    std::fill(
        dist.begin(), dist.end(),
        std::make_pair(std::numeric_limits<int64_t>::max(), NetworKit::index{0}));
    std::fill(visited.begin(), visited.end(), false);
    dist[source] = {0, 0};
    pq.push({0, source});

    while (!pq.empty()) {
        auto [d, u] = pq.top();
        pq.pop();

        if (visited[u]) continue;
        visited[u] = true;

        for (NetworKit::index edge_idx : neighbors[u]) {
            const Edge& edge = edges[edge_idx];
            if (edge.capacity >= edge.flow + delta) {
                node v = edge.to;
                int64_t new_dist = d + edge.cost - potential[u] + potential[v];
                if (new_dist < dist[v].first) {
                    dist[v] = {new_dist, edge_idx};
                    if (is_added_uncapacitated(v)) {
                        for (NetworKit::index e2 : neighbors[v]) {
                            const Edge& edge2 = edges[e2];
                            if (edge2.to != u && edge2.capacity >= edge2.flow + delta) {
                                node k = edge2.to;
                                int64_t nd2 = new_dist + edge2.cost - potential[v] + potential[k];
                                if (nd2 < dist[k].first) {
                                    dist[k] = {nd2, e2};
                                    pq.push({nd2, k});
                                }
                            }
                        }
                    } else {
                        pq.push({new_dist, v});
                    }
                }
            }
        }
    }
}

void OrlinMinimumCostFlow::make_reduced_costs_nonnegative() {
    constexpr int64_t INF = std::numeric_limits<int32_t>::max();
    std::vector<int64_t> d(max_nodeid, INF);
    d[0] = 0;

    bool changed = true;
    for (NetworKit::count iter = 0; iter < nodes_number && changed; ++iter) {
        changed = false;
        for (const Edge& edge : edges) {
            if (edge.from == edge.to || edge.capacity <= edge.flow) continue;
            if (d[edge.from] == INF) continue;
            int64_t len = edge.cost - potential[edge.from] + potential[edge.to];
            if (d[edge.from] + len < d[edge.to]) {
                d[edge.to] = d[edge.from] + len;
                changed = true;
            }
        }
    }

    for (node v = 0; v < max_nodeid; ++v) {
        if (d[v] == INF) continue;
        potential[v] -= d[v];
        potential_computed[v] -= d[v];
    }
}

void OrlinMinimumCostFlow::compute_final_flows() {
    NetworKit::Graph maxflow_graph(max_nodeid, true, true);
    original_graph.forEdges([&](node u, node v) {
        auto cost = network.cost[{u, v}];

        if (cost - potential_computed[u] + potential_computed[v] == 0) {
            NetworKit::edgeweight max = std::numeric_limits<int>::max();
            maxflow_graph.addEdge(u, v, max);
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
