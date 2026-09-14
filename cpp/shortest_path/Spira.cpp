#include "shortest_path/Spira.hpp"

#include <limits>
#include <queue>
#include <stdexcept>

#include <networkit/graph/GraphTools.hpp>

namespace Koala {

namespace {

static const NetworKit::edgeweight INF = std::numeric_limits<NetworKit::edgeweight>::max();

struct Edge {
    NetworKit::node target;
    NetworKit::edgeweight weight;
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

void SpiraAPSP::checkInput() const {
    if (!graph->isWeighted()) {
        throw std::invalid_argument("Graph must be weighted for Spira's algorithm.");
    }
}

void SpiraAPSP::run() {
    const NetworKit::count n = graph->upperNodeIdBound();
    const NetworKit::count actual_nodes = graph->numberOfNodes();

    NetworKit::GraphTools::sortEdgesByWeight(*graph);

    std::vector<std::vector<Edge>> adj(n);
    for (const auto u : graph->nodeRange()) {
        for (const auto [v, w] : graph->weightNeighborRange(u)) {
            adj[u].push_back({v, w});
        }
    }

    distances.assign(n, std::vector<NetworKit::edgeweight>(n, INF));

    for (const auto origin : graph->nodeRange()) {
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
    }

    computeDiameter();

    hasRun = true;
}

}  // namespace Koala
