#include <vector>
#include <cstdlib>
#include <ctime>
#include <climits>
#include <random>

#include <min_cut/KargerMinCut.hpp>

namespace Koala {

int KargerMinCut::find(std::vector<int>& parent, int i) {
    if (parent[i] == i) return i;

    return find(parent, parent[i]);
}

void KargerMinCut::unionSub(std::vector<int>& parent, std::vector<int>& rank, int x, int y) {
    int xRoot = find(parent, x);
    int yRoot = find(parent, y);

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
    int vertices = graph->numberOfNodes();
    std::vector<int> parent(vertices), rank(vertices, 0);
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

    for (int i = 0; i < graph->numberOfNodes(); ++i) {
        parent[i] = i;
    }

    while (vertices > 2) {
        double pick = dis(gen);
        double cumulativeWeight = 0;

        int selectedEdge = -1;
        for (int i = 0; i < edges.size(); ++i) {
            cumulativeWeight += edges[i].weight;
            if (cumulativeWeight >= pick) {
                selectedEdge = i;
                break;
            }
        }

        int subset1 = find(parent, edges[selectedEdge].u);
        int subset2 = find(parent, edges[selectedEdge].v);

        if (subset1 != subset2) {
            unionSub(parent, rank, subset1, subset2);
            vertices--;

            totalWeight -= edges[selectedEdge].weight;
        }
    }

    graph->forEdges([&](NetworKit::node u, NetworKit::node v, double weight) {
        int subset1 = find(parent, u);
        int subset2 = find(parent, v);
        if (subset1 != subset2) {
            minCutValue += weight;
        }
    });
}


}  // namespace Koala
