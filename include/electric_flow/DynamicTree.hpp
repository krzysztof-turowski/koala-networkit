#pragma once

#include <tuple>
#include <vector>
#include <map>

#include "networkit/graph/Graph.hpp"

namespace Koala {

class DynamicTree {
public:
    DynamicTree(int n, std::vector<std::vector<double>>& weights);

    void link(int u, int v);
    void cut(int u, int v);
    int findRoot(int v);
    void pathAdd(int u, int v, double c);
    std::pair<int, int> pathMin(int u, int v);
    double pathSum(int u, int v);

    std::vector<std::vector<double>>& weights;
    NetworKit::Graph graph;
};

}; // namespace koala
