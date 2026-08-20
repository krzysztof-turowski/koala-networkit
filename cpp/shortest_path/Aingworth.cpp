#include "shortest_path/Aingworth.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <stdexcept>

namespace Koala {

static const int INF = std::numeric_limits<int>::max();

AingworthAlgorithm::AingworthAlgorithm(const NetworKit::Graph& G)
    : G(&G), diameter(0) {
    if (G.isDirected()) {
        throw std::invalid_argument("Graph must be undirected for Aingworth's algorithm.");
    }
    if (G.isWeighted()) {
        throw std::invalid_argument("Graph must be unweighted for Aingworth's algorithm.");
    }
}

std::vector<int> AingworthAlgorithm::bfs(NetworKit::node src, const std::vector<std::vector<NetworKit::node>>& adj, NetworKit::count n) {
    std::vector<int> dist(n, INF);
    dist[src] = 0;
    std::queue<NetworKit::node> q;
    q.push(src);
    while (!q.empty()) {
        NetworKit::node u = q.front();
        q.pop();
        for (NetworKit::node v : adj[u]) {
            if (dist[v] == INF) {
                dist[v] = dist[u] + 1;
                q.push(v);
            }
        }
    }
    return dist;
}

std::vector<int> AingworthAlgorithm::bfsLow(NetworKit::node src, const std::vector<std::vector<NetworKit::node>>& adj, const std::vector<bool>& is_low, NetworKit::count n) {
    std::vector<int> dist(n, INF);
    dist[src] = 0;
    std::queue<NetworKit::node> q;
    q.push(src);
    while (!q.empty()) {
        NetworKit::node u = q.front();
        q.pop();
        for (NetworKit::node v : adj[u]) {
            if (is_low[v] && dist[v] == INF) {
                dist[v] = dist[u] + 1;
                q.push(v);
            }
        }
    }
    return dist;
}

std::vector<NetworKit::node> AingworthAlgorithm::dominatingSet(const std::vector<NetworKit::node>& high, const std::vector<std::vector<NetworKit::node>>& adj, const std::vector<bool>& is_high, NetworKit::count n) {
    std::vector<bool> dominated(n, false);
    std::vector<NetworKit::node> D;
    for (NetworKit::node v : high) {
        if (!dominated[v]) {
            D.push_back(v);
            dominated[v] = true;
            for (NetworKit::node u : adj[v]) {
                if (is_high[u]) {
                    dominated[u] = true;
                }
            }
        }
    }
    return D;
}

void AingworthAlgorithm::run() {
    const NetworKit::count n = G->upperNodeIdBound();
    const double nn = static_cast<double>(G->numberOfNodes());
    const int s = std::max(1, static_cast<int>(std::ceil(std::sqrt(nn * std::log2(nn)))));

    std::vector<std::vector<NetworKit::node>> adj(n);
    G->forNodes([&](NetworKit::node u) {
        G->forNeighborsOf(u, [&](NetworKit::node v) {
            adj[u].push_back(v);
        });
    });

    std::vector<bool> is_low(n, false);
    std::vector<bool> is_high(n, false);
    std::vector<NetworKit::node> L, H;

    G->forNodes([&](NetworKit::node v) {
        if (static_cast<int>(adj[v].size()) < s) {
            is_low[v] = true;
            L.push_back(v);
        } else {
            is_high[v] = true;
            H.push_back(v);
        }
    });

    distances.assign(n, std::vector<int>(n, INF));
    G->forNodes([&](NetworKit::node v) {
        distances[v][v] = 0;
    });

    auto D = dominatingSet(H, adj, is_high, n);
    std::vector<bool> in_D(n, false);
    for (NetworKit::node v : D) {
        in_D[v] = true;
    }

    for (NetworKit::node v : D) {
        auto dist = bfs(v, adj, n);
        for (NetworKit::count u = 0; u < n; ++u) {
            if (dist[u] < distances[v][u]) {
                distances[v][u] = dist[u];
                distances[u][v] = dist[u];
            }
        }
    }

    for (NetworKit::node v : L) {
        auto dist = bfsLow(v, adj, is_low, n);
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

    diameter = 0;
    G->forNodes([&](NetworKit::node u) {
        G->forNodes([&](NetworKit::node v) {
            if (distances[u][v] < INF) {
                diameter = std::max(diameter, distances[u][v]);
            }
        });
    });

    hasRun = true;
}

const std::vector<std::vector<int>>& AingworthAlgorithm::getDistances() const {
    assureFinished();
    return distances;
}

int AingworthAlgorithm::getDistance(NetworKit::node u, NetworKit::node v) const {
    assureFinished();
    return distances[u][v];
}

int AingworthAlgorithm::getDiameter() const {
    assureFinished();
    return diameter;
}

}
