#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <algorithm>
#include <climits>

#include <min_cut/HaoOrlinMinCut.hpp>
#include <min_cut/PushRelabelMinCut.hpp>

namespace Koala {

void HaoOrlinMinCut::findMinimumCut() {
    bestValue = INT_MAX;
    cut.clear();
    std::vector<bool> visited(numVertices, false);
    visited[source] = true;
    std::vector<int> S = {source};

    while (S.size() < numVertices) {
        int t_prime = -1;
        for (int i = 0; i < numVertices; ++i) {
            if (!visited[i]) {
                t_prime = i;
                break;
            }
        }

        PushRelabelMinCut minCutSolver(numVertices, source, t_prime, graph);
        minCutSolver.solve();
        std::vector<int> currentCut = minCutSolver.getMinCut();
        int z = minCutSolver.getMaxFlow();

        if (z < bestValue) {
            bestValue = z;
            cut = currentCut;
        }

        visited[t_prime] = true;
        S.push_back(t_prime);
    }
}

HaoOrlinMinCut::HaoOrlinMinCut(int vertices, int src, int snk, const std::vector<std::vector<int>>& graphMatrix)
    : numVertices(vertices), source(src), sink(snk), graph(graphMatrix) {
}


long long HaoOrlinMinCut::getMinCut() const {
    return bestValue;
}

}  // namespace Koala
