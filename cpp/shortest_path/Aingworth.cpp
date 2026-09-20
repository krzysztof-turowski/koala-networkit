#include "shortest_path/Aingworth.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <stdexcept>
#include <utility>

#include <networkit/distance/BFS.hpp>

namespace Koala {

static const int INF = std::numeric_limits<int>::max();

void AingworthAPSP::checkInput() const {
    if (graph->isDirected()) {
        throw std::invalid_argument("Graph must be undirected for Aingworth's algorithm.");
    }
    if (graph->isWeighted()) {
        throw std::invalid_argument("Graph must be unweighted for Aingworth's algorithm.");
    }
}

std::vector<int> AingworthAPSP::bfs(NetworKit::node src, const NetworKit::Graph& graph, NetworKit::count n) {
    NetworKit::BFS bfs(graph, src, false, false);
    bfs.run();
    const auto& d = bfs.getDistances();

    std::vector<int> dist(n, INF);
    for (NetworKit::node u = 0; u < n; ++u) {
        if (d[u] < static_cast<double>(INF)) {
            dist[u] = static_cast<int>(d[u]);
        }
    }
    return dist;
}

std::vector<NetworKit::node> AingworthAPSP::dominatingSet(const std::vector<std::vector<NetworKit::node>>& adj, const std::vector<bool>& is_high, NetworKit::count n) {
    std::vector<int> cover_count(n, 0);
    for (NetworKit::node v = 0; v < n; ++v) {
        int c = is_high[v] ? 1 : 0;
        for (NetworKit::node u : adj[v]) {
            if (is_high[u]) {
                ++c;
            }
        }
        cover_count[v] = c;
    }

    std::vector<bool> covered(n, false);
    std::vector<bool> in_D(n, false);
    std::vector<NetworKit::node> D;

    using Item = std::pair<int, NetworKit::node>;
    std::priority_queue<Item> pq;
    for (NetworKit::node v = 0; v < n; ++v) {
        if (cover_count[v] > 0) {
            pq.push({cover_count[v], v});
        }
    }

    auto markCovered = [&](NetworKit::node w) {
        covered[w] = true;
        --cover_count[w];
        if (cover_count[w] > 0) {
            pq.push({cover_count[w], w});
        }
        for (NetworKit::node x : adj[w]) {
            --cover_count[x];
            if (cover_count[x] > 0) {
                pq.push({cover_count[x], x});
            }
        }
    };

    while (!pq.empty()) {
        auto [cnt, v] = pq.top();
        pq.pop();
        if (in_D[v] || cnt != cover_count[v]) {
            continue;  // already chosen, or a stale heap entry
        }
        if (cnt <= 0) {
            break;  // nothing uncovered remains
        }
        in_D[v] = true;
        D.push_back(v);
        if (is_high[v] && !covered[v]) {
            markCovered(v);
        }
        for (NetworKit::node u : adj[v]) {
            if (is_high[u] && !covered[u]) {
                markCovered(u);
            }
        }
    }
    return D;
}

void AingworthAPSP::run() {
    const NetworKit::count n = graph->upperNodeIdBound();
    const double nn = static_cast<double>(graph->numberOfNodes());
    const int s = threshold > 0
        ? threshold
        : std::max(1, static_cast<int>(std::ceil(std::sqrt(nn * std::log2(nn)))));

    std::vector<std::vector<NetworKit::node>> adj(n);
    for (const auto u : graph->nodeRange()) {
        for (const auto v : graph->neighborRange(u)) {
            adj[u].push_back(v);
        }
    }

    std::vector<bool> is_low(n, false);
    std::vector<bool> is_high(n, false);
    std::vector<NetworKit::node> L;

    for (const auto v : graph->nodeRange()) {
        if (static_cast<int>(graph->degree(v)) < s) {
            is_low[v] = true;
            L.push_back(v);
        } else {
            is_high[v] = true;
        }
    }

    NetworKit::Graph lowG(n, false, false);
    for (const auto v : graph->nodeRange()) {
        if (is_high[v]) {
            lowG.removeNode(v);
        }
    }

    for (NetworKit::node u : L) {
        for (NetworKit::node v : adj[u]) {
            if (is_low[v] && u < v) {
                lowG.addEdge(u, v);
            }
        }
    }

    distances.assign(n, std::vector<int>(n, INF));
    for (const auto v : graph->nodeRange()) {
        distances[v][v] = 0;
    }

    auto D = dominatingSet(adj, is_high, n);
    std::vector<bool> in_D(n, false);
    for (NetworKit::node v : D) {
        in_D[v] = true;
    }

    for (NetworKit::node v : D) {
        auto dist = bfs(v, *graph, n);
        for (NetworKit::count u = 0; u < n; ++u) {
            if (dist[u] < distances[v][u]) {
                distances[v][u] = dist[u];
                distances[u][v] = dist[u];
            }
        }
    }

    for (NetworKit::node v : L) {
        auto dist = bfs(v, lowG, n);
        for (NetworKit::count u = 0; u < n; ++u) {
            if (dist[u] < INF && dist[u] < distances[v][u]) {
                distances[v][u] = dist[u];
                distances[u][v] = dist[u];
            }
        }
    }

    for (NetworKit::count u = 0; u < n; ++u) {
        if (in_D[u]) continue;
        for (NetworKit::count v = u + 1; v < n; ++v) {
            if (in_D[v]) continue;
            for (NetworKit::node w : D) {
                if (distances[w][u] < INF && distances[w][v] < INF) {
                    int via_w = distances[w][u] + distances[w][v];
                    if (via_w < distances[u][v]) {
                        distances[u][v] = via_w;
                        distances[v][u] = via_w;
                    }
                }
            }
        }
    }

    computeDiameter();

    hasRun = true;
}

}  // namespace Koala
