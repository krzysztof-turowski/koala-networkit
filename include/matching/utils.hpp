#pragma once

#include <cmath>
#include <random>
#include <iostream>
#include <map>

#include <networkit/graph/Graph.hpp>

using namespace std;
using namespace NetworKit;

namespace Koala {
    constexpr double EPS = 1e-8;

    inline bool eq(double a, double b) {
        return fabs(a - b) <= EPS;
    }

    inline double generateRandom() {
        return (double)rand() / RAND_MAX; //TODO
    }

    inline void printGraph(const Graph& G) {
        cout << "v\t";
        G.forNodes([](node v) {
            cout << v << ' ';
            });
        cout << '\n';
        G.forEdges([](node u, node v) {
            cout << "e\t" << u << ' ' << v << endl;
            });
    }

    inline pair<Graph, std::vector<int>> reindexGraph(const Graph& G) {
        int n = G.numberOfNodes();
        Graph G1(n, false, false);
        map<int, int> indexes;
        vector<int> labels(n);
        for (int i=0; auto v: G.nodeRange()) {
            indexes[v] = i;
            labels[i] = v;
            i++;
        }
        G.forEdges([&](node u, node v) { G1.addEdge(indexes[u], indexes[v]); });
        return { G1, labels };
    }
}