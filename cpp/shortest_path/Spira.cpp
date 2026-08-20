#include "shortest_path/Spira.hpp"

#include <algorithm>
#include <limits>
#include <queue>

namespace Koala {

namespace {

static const NetworKit::edgeweight INF = std::numeric_limits<NetworKit::edgeweight>::max();

struct Edge {
    NetworKit::node target;
    NetworKit::edgeweight weight;

    bool operator<(const Edge& other) const {
        return weight < other.weight;
    }
};

struct State {
    NetworKit::edgeweight dist;
    NetworKit::node from;
    NetworKit::index edge_idx;

    bool operator>(const State& other) const {
        return dist > other.dist;
    }
};

}

SpiraAlgorithm::SpiraAlgorithm(const NetworKit::Graph& G) : G(&G) {
    if (!G.isWeighted()) {
        throw std::invalid_argument("Graph must be weighted for Spira's algorithm.");
    }
}

void SpiraAlgorithm::run() {
    const NetworKit::count n = G->upperNodeIdBound();
    const NetworKit::count actual_nodes = G->numberOfNodes();

    std::vector<std::vector<Edge>> adj(n);
    G->forNodes([&](NetworKit::node u) {
        G->forNeighborsOf(u, [&](NetworKit::node v, NetworKit::edgeweight w) {
            adj[u].push_back({v, w});
        });
        std::sort(adj[u].begin(), adj[u].end());
    });

    distances.assign(n, std::vector<NetworKit::edgeweight>(n, INF));

    G->forNodes([&](NetworKit::node origin) {
        distances[origin][origin] = 0.0;
        std::vector<bool> labeled(n, false);
        labeled[origin] = true;
        NetworKit::count count_labeled = 1;

        std::priority_queue<State, std::vector<State>, std::greater<State>> pq;
        NetworKit::node current = origin;

        while (count_labeled < actual_nodes) {
            if (!adj[current].empty()) {
                pq.push({distances[origin][current] + adj[current][0].weight, current, 0});
            }

            bool found = false;
            while (!pq.empty() && !found) {
                auto [d, u, idx] = pq.top();
                pq.pop();

                if (idx + 1 < static_cast<NetworKit::index>(adj[u].size())) {
                    pq.push({distances[origin][u] + adj[u][idx + 1].weight, u, idx + 1});
                }

                NetworKit::node v = adj[u][idx].target;
                if (!labeled[v]) {
                    labeled[v] = true;
                    distances[origin][v] = d;
                    current = v;
                    count_labeled++;
                    found = true;
                }
            }

            if (!found) break;
        }
    });

    hasRun = true;
}

const std::vector<std::vector<NetworKit::edgeweight>>& SpiraAlgorithm::getDistances() const {
    assureFinished();
    return distances;
}

NetworKit::edgeweight SpiraAlgorithm::getDistance(NetworKit::node u, NetworKit::node v) const {
    assureFinished();
    return distances[u][v];
}

}
