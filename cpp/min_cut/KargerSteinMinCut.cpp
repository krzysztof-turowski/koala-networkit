#include <cstdlib>
#include <ctime>
#include <cmath>
#include <iostream>
#include <climits>

#include <min_cut/KargerSteinMinCut.hpp>

namespace Koala {

KargerSteinMinCut::KargerSteinMinCut(int vertices, int repeats, const std::vector<std::vector<int>>& graphMatrix) : numVertices(vertices), numRepeats(repeats) {
    srand(time(NULL));
    
    for (int i = 0; i < vertices; ++i) {
        for (int j = 0; j < vertices; ++j) {
            if (graphMatrix[i][j] > 0) {
                addEdge(i, j, graphMatrix[i][j]);
            }
        }
    }
}

void KargerSteinMinCut::addEdge(int u, int v, int weight) {
    for (int i = 0; i < weight; ++i) { // Adding multiple edges for the given weight
        edges.push_back({u, v});
    }
}

int KargerSteinMinCut::find(std::vector<int>& parent, int i) {
    if (parent[i] == i)
        return i;
    return find(parent, parent[i]);
}

void KargerSteinMinCut::unionSub(std::vector<int>& parent, std::vector<int>& rank, int x, int y) {
    int xRoot = find(parent, x);
    int yRoot = find(parent, y);

    if (rank[xRoot] < rank[yRoot])
    {
        parent[xRoot] = yRoot;
    }
    else if (rank[xRoot] > rank[yRoot])
    {
        parent[yRoot] = xRoot;
    }
    else
    {
        parent[yRoot] = xRoot;
        rank[xRoot]++;
    }
}

int KargerSteinMinCut::recursiveMinCut(int currentV) {
    if (currentV <= 6)
    {
        return standardMinCut();
    }

    int T = (int)std::ceil(currentV / std::sqrt(2));
    std::vector<int> parent(currentV), rank(currentV, 0);
    for (int i = 0; i < currentV; ++i) {
        parent[i] = i;
    }

    int vertices = currentV;
    while (vertices > T) {
        int i = rand() % edges.size();
        int subset1 = find(parent, edges[i].u);
        int subset2 = find(parent, edges[i].v);

        if (subset1 != subset2) {
            unionSub(parent, rank, subset1, subset2);
            vertices--;
        }
    }

    return std::min(recursiveMinCut(T), recursiveMinCut(T));
}

void KargerSteinMinCut::solve() {
    bestMinCut = INT_MAX;
    for (int i = 0; i < numRepeats; ++i) {
        int currentMinCut = recursiveMinCut(numVertices);

        bestMinCut = std::min(bestMinCut, currentMinCut);
    }
}

int KargerSteinMinCut::getMinCut() {
    return bestMinCut;
}

int KargerSteinMinCut::standardMinCut() {
    std::vector<int> parent(numVertices), rank(numVertices, 0);
    int vertices = numVertices;
    int minCut = 0;
    
    for (int i = 0; i < numVertices; ++i) {
        parent[i] = i;
    }

    while (vertices > 2) {
        int i = rand() % edges.size();
        int subset1 = find(parent, edges[i].u);
        int subset2 = find(parent, edges[i].v);

        if (subset1 != subset2) {
            unionSub(parent, rank, subset1, subset2);
            vertices--;
        }
    }

    for (auto& edge : edges) {
        int subset1 = find(parent, edge.u);
        int subset2 = find(parent, edge.v);
        if (subset1 != subset2)
            minCut++;
    }
    return minCut / 2; // Each edge counted twice
}

} // namespace Koala
