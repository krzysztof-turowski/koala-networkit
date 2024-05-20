#include <cstdlib>
#include <ctime>
#include <cmath>
#include <iostream>
#include <climits>

#include <min_cut/KargerSteinMinCut.hpp>

namespace Koala {

int KargerSteinMinCut::find(std::vector<int>& parent, int i) {
    if (parent[i] == i)
        return i;
    return find(parent, parent[i]);
}

void KargerSteinMinCut::unionSub(std::vector<int>& parent, std::vector<int>& rank, int x, int y) {
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

int KargerSteinMinCut::recursiveMinCut(int currentV) {
    if (currentV <= 6) {
        return standardMinCut();
    }

    int T = static_cast<int>(std::ceil(currentV / std::sqrt(2)));
    std::vector<int> parent(currentV), rank(currentV, 0);
    for (int i = 0; i < currentV; ++i) {
        parent[i] = i;
    }

    std::random_device rd;
    std::mt19937 gen(rd());

    double totalWeight = 0;
    std::vector<NetworKit::WeightedEdge> weightedEdges;
    graph->forEdges([&](NetworKit::node u, NetworKit::node v, double weight) {
        weightedEdges.emplace_back(u, v, weight);
        totalWeight += weight;
    });

    int vertices = currentV;
    std::uniform_real_distribution<> dis(0.0, totalWeight);

    while (vertices > T) {
        double pick = dis(gen);
        double cumulativeWeight = 0;

        int selectedEdge = -1;
        for (int i = 0; i < weightedEdges.size(); ++i) {
            cumulativeWeight += weightedEdges[i].weight;
            if (cumulativeWeight >= pick) {
                selectedEdge = i;
                break;
            }
        }

        int subset1 = find(parent, weightedEdges[selectedEdge].u);
        int subset2 = find(parent, weightedEdges[selectedEdge].v);

        if (subset1 != subset2) {
            unionSub(parent, rank, subset1, subset2);
            vertices--;

            totalWeight -= weightedEdges[selectedEdge].weight;
        }
    }

    return std::min(recursiveMinCut(T), recursiveMinCut(T));
}

int KargerSteinMinCut::standardMinCut() {
    int vertices = graph->numberOfNodes();
    std::vector<int> parent(vertices), rank(vertices, 0);
    double minCut = 0;

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
            minCut += weight;
        }
    });

    return minCut;
}

void KargerSteinMinCut::run() {
    double bestMinCutValue = INT_MAX;
    for (int i = 0; i < repeat; ++i) {
        runOnce();

        bestMinCutValue = std::min(bestMinCutValue, minCutValue);
    }
    minCutValue = bestMinCutValue;
}

void KargerSteinMinCut::runOnce() {
    minCutValue = recursiveMinCut(graph->numberOfNodes());
}

}  // namespace Koala
