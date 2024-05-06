#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <climits>

#include <min_cut/KargerMinCut.hpp>

namespace Koala {

KargerMinCut::KargerMinCut(int vertices, int repeat, const std::vector<std::vector<int>>& graphMatrix) : numVertices(vertices), numRepeat(repeat) {
    srand(time(NULL));
    
    for (int i = 0; i < vertices; ++i) {
        for (int j = 0; j < vertices; ++j) {
            if (graphMatrix[i][j] > 0) {
                addEdge(i, j, graphMatrix[i][j]);
            }
        }
    }
}

void KargerMinCut::addEdge(int u, int v, int weight) {
    for (int i = 0; i < weight; ++i) { // Adding multiple edges for the given weight
        edges.push_back({u, v});
    }
}

int KargerMinCut::find(std::vector<int>& parent, int i) {
    if (parent[i] == i) return i;

    return find(parent, parent[i]);
}

void KargerMinCut::unionSub(std::vector<int>& parent, std::vector<int>& rank, int x, int y) {
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


void KargerMinCut::solve() {
    bestMinCut = INT_MAX;
    for (int i = 0; i < numRepeat; ++i) {
        int currentMinCut = findMinCut();

        bestMinCut = std::min(bestMinCut, currentMinCut);
    }
}

int KargerMinCut::findMinCut() {
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

int KargerMinCut::getMinCut() {
    return bestMinCut;
}

} // namespace Koala