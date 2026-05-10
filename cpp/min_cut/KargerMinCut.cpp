#include <algorithm>
#include <climits>
#include <cstdlib>
#include <ctime>
#include <random>
#include <vector>

#include <min_cut/KargerMinCut.hpp>

namespace Koala {

NetworKit::node KargerMinCut::find(std::vector<NetworKit::node>& parent, NetworKit::node i) {
    if (parent[i] == i) return i;

    return find(parent, parent[i]);
}

void KargerMinCut::unionSub(
        std::vector<NetworKit::node>& parent, std::vector<NetworKit::count>& rank,
        NetworKit::node x, NetworKit::node y) {
    auto xRoot = find(parent, x);
    auto yRoot = find(parent, y);

    if (rank[xRoot] < rank[yRoot]) {
        parent[xRoot] = yRoot;
    } else if (rank[xRoot] > rank[yRoot]) {
        parent[yRoot] = xRoot;
    } else {
        parent[yRoot] = xRoot;
        rank[xRoot]++;
    }
}


void KargerMinCut::run() {
    double bestMinCutValue = INT_MAX;
    for (int i = 0; i < repeat; ++i) {
        runOnce();

        bestMinCutValue = std::min(bestMinCutValue, minCutValue);
    }
    minCutValue = bestMinCutValue;
}

void KargerMinCut::runOnce() {
    auto vertices = graph->numberOfNodes();
    std::vector<NetworKit::node> parent(vertices);
    std::vector<NetworKit::count> rank(vertices, 0);
    minCutValue = 0;

    std::random_device rd;
    std::mt19937 gen(rd());

    double totalWeight = 0;
    std::vector<NetworKit::WeightedEdge> edges;
    graph->forEdges([&](NetworKit::node u, NetworKit::node v, double weight) {
        edges.emplace_back(u, v, weight);
        totalWeight += weight;
    });

    std::uniform_real_distribution<> dis(0.0, totalWeight);

    for (auto node : graph->nodeRange()) {
        parent[node] = node;
    }

    while (vertices > 2) {
        double pick = dis(gen);
        double cumulativeWeight = 0;

        std::size_t selectedEdge = -1;
        for (std::size_t i = 0; i < edges.size(); ++i) {
            cumulativeWeight += edges[i].weight;
            if (cumulativeWeight >= pick) {
                selectedEdge = i;
                break;
            }
        }

        auto subset1 = find(parent, edges[selectedEdge].u);
        auto subset2 = find(parent, edges[selectedEdge].v);

        if (subset1 != subset2) {
            unionSub(parent, rank, subset1, subset2);
            vertices--;

            totalWeight -= edges[selectedEdge].weight;
        }
    }

    graph->forEdges([&](NetworKit::node u, NetworKit::node v, double weight) {
        auto subset1 = find(parent, u);
        auto subset2 = find(parent, v);
        if (subset1 != subset2) {
            minCutValue += weight;
        }
    });
}

}  // namespace Koala
