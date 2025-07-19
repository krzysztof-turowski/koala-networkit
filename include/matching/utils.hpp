#pragma once

#include <cmath>
#include <random>
#include <iostream>
#include <map>

#include <networkit/graph/Graph.hpp>
#include <Eigen/Core>

using namespace std;
using namespace NetworKit;
using namespace Eigen;

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

    inline tuple<Graph, std::vector<int>, MatrixXd> reindexGraph(const Graph& G, const MatrixXd& AG) {
        int n = G.numberOfNodes();
        auto [G1, labels] = reindexGraph(G);

        Eigen::MatrixXd AG1(n, n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                AG1(i, j) = AG(labels[i], labels[j]);
            }
        }
        
        return { G1, labels, AG1 };
    }
}