#include <iostream>
#include <vector>
#include <algorithm>
#include <climits>
#include <queue>

#include <min_cut/StoerWagnerMinCut.hpp>

namespace Koala {

StoerWagnerMinCut::StoerWagnerMinCut(int vertices, const std::vector<std::vector<int>>& graphMatrix)
    : numVertices(vertices), graph(graphMatrix), vertices(vertices), bestMinCut(INT_MAX) {
    for (int i = 0; i < vertices; ++i) {
        this->vertices[i] = i;
    }
}

void StoerWagnerMinCut::solve() {
    bestMinCut = findMinCut(graph);
}

int StoerWagnerMinCut::getMinCut() const {
    return bestMinCut;
}

int StoerWagnerMinCut::findMinCut(std::vector<std::vector<int>>& weights) {
    std::vector<int> localVertices = vertices;
    int minCut = INT_MAX;
    while (localVertices.size() > 1) {
        minCut = std::min(minCut, minCutPhase(localVertices, weights));
    }
    return minCut;
}

int StoerWagnerMinCut::minCutPhase(std::vector<int>& vertices, std::vector<std::vector<int>>& weights) {
    int size = vertices.size();
    std::vector<int> weightsPhase(size, 0);
    std::vector<bool> added(size, false);

    int prev, last = 0;

    for (int i = 0; i < size; ++i) {
        prev = last;
        last = -1;
        for (int j = 0; j < size; ++j) {
            if (!added[vertices[j]] && (last == -1 || weightsPhase[vertices[j]] > weightsPhase[vertices[last]])) {
                last = j;
            }
        }

        if (i == size - 1) {
            int minCut = weightsPhase[vertices[last]];
            for (int j = 0; j < size; ++j) {
                weights[vertices[prev]][vertices[j]] += weights[vertices[last]][vertices[j]];
                weights[vertices[j]][vertices[prev]] += weights[vertices[j]][vertices[last]];
            }
            vertices.erase(vertices.begin() + last);
            return minCut;
        }
        added[vertices[last]] = true;
        for (int j = 0; j < size; ++j) {
            if (!added[vertices[j]]) {
                weightsPhase[vertices[j]] += weights[vertices[last]][vertices[j]];
            }
        }
    }
    return INT_MAX;
}

} // namespace Koala
