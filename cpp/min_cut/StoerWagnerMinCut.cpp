#include <iostream>
#include <vector>
#include <algorithm>
#include <climits>
#include <queue>

#include <min_cut/StoerWagnerMinCut.hpp>

namespace Koala {

void StoerWagnerMinCut::run() {
    std::vector<NetworKit::node> localVertices(graph->numberOfNodes());
    std::vector<std::vector<double>> edges(graph->numberOfNodes(),
            std::vector<double>(graph->numberOfNodes(), 0));
    minCutValue = INT_MAX;

    graph->forNodes([&](NetworKit::node u) {
        localVertices[u] = u;
    });

    graph->forEdges([&](NetworKit::node u, NetworKit::node v, double w) {
        edges[u][v] = w;
        edges[v][u] = w;
    });

    while (localVertices.size() > 1) {
        minCutValue = std::min(minCutValue, minCutPhase(localVertices, edges));
    }
}

double StoerWagnerMinCut::minCutPhase(std::vector<NetworKit::node>& vertices,
        std::vector<std::vector<double>>& edges) {
    int size = vertices.size();
    std::vector<double> weightsPhase(size, 0);
    std::vector<bool> added(size, false);

    int prev, last = 0;

    for (int i = 0; i < size; ++i) {
        prev = last;
        last = -1;
        for (int j = 0; j < size; ++j) {
            if (!added[j] && (last == -1 || weightsPhase[j] > weightsPhase[last])) {
                last = j;
            }
        }

        if (i == size - 1) {
            double minCut = weightsPhase[last];
            for (int j = 0; j < size; ++j) {
                if (j != last) {
                    edges[vertices[prev]][vertices[j]] += edges[vertices[last]][vertices[j]];
                    edges[vertices[j]][vertices[prev]] += edges[vertices[j]][vertices[last]];
                }
            }
            vertices.erase(vertices.begin() + last);
            return minCut;
        }
        added[last] = true;
        for (int j = 0; j < size; ++j) {
            if (!added[j]) {
                weightsPhase[j] += edges[vertices[last]][vertices[j]];
            }
        }
    }
    return INT_MAX;
}

}  // namespace Koala
