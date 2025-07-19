#include "electric_flow/DynamicTree.hpp"

#include <algorithm>
#include <stack>
#include <cassert>

using namespace std;
using namespace NetworKit;

namespace Koala {

vector<int> getPath(const Graph& graph, int u, int v);

DynamicTree::DynamicTree(int n, vector<vector<double>>& weights): graph(Graph(n)), weights(weights) {}

void DynamicTree::link(int u, int v) {
    assert(findRoot(u) != findRoot(v));
    graph.addEdge(u, v);
}

void DynamicTree::cut(int u, int v) {
    if (graph.hasEdge(u,v)) {
    graph.removeEdge(u, v);
    }
}

int DynamicTree::findRoot(int v) {
    int root = v;
    vector<int> s;
    vector<bool> visited(graph.numberOfNodes(), false);
    s.push_back(v);
    while (!s.empty()) {
        int x = s.back();
        s.pop_back();
        visited[x] = true;
        root = min(root, x);
        for (auto y: graph.neighborRange(x)) {
            if (!visited[y]) {
                s.push_back(y);
            }
        }
    }
    return root;
}

void DynamicTree::pathAdd(int u, int v, double c) {
    auto path = getPath(graph, u, v);
    for (int i = 1; i < path.size(); ++i) {
        weights[path[i-1]][path[i]] -= c;
        weights[path[i]][path[i-1]] += c;
    }
}

std::pair<int, int> DynamicTree::pathMin(int u, int v) {
    auto path = getPath(graph, u, v);
    double m = 1e+37;
    std::pair<int, int> mv = {-1,-1};
    for (int i = 1; i < path.size(); ++i) {
        double w = weights[path[i-1]][path[i]];
        if (w >= 0 && w <= m) {
            m = w;
            mv = {path[i-1], path[i]};
        }
    }
    return mv;
}

double DynamicTree::pathSum(int u, int v) {
    auto path = getPath(graph, u, v);
    double s = 0;
    for (int i = 1; i < path.size(); ++i) {
        s += weights[path[i-1]][path[i]];
    }
    return s;
}


vector<int> getPath(const Graph& graph, int u, int v) {
    vector<int> s;
    vector<int> parent(graph.numberOfNodes(), -1);
    s.push_back(u);
    parent[u] = u;
    while (!s.empty()) {
        int x = s.back();
        if (x == v) {
            break;
        }
        s.pop_back();
        for (auto y: graph.neighborRange(x)) {
            if (parent[y] == -1) {
                parent[y] = x;
                s.push_back(y);
            }
        }
    }
    vector<int> path;
    int x = v;
    while (x != u) {
        path.push_back(x);
        x = parent[x];
    }
    path.push_back(u);

    reverse(path.begin(), path.end());

    return path;
}

} // namespace Koala