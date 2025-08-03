#include <networkit/graph/Graph.hpp>

#include <matching/gaussian_matching/DynamicComponents.hpp>
#include <matching/gaussian_matching/utils.hpp>

using namespace NetworKit;
using namespace std;

namespace Koala {
    void dfs(const Graph& G, vector<bool>& visited, int v);

    DynamicComponents::DynamicComponents(const Graph& G) : G(G) {}

    void DynamicComponents::addEdge(int u, int v) {
        G.addEdge(u, v);
    }

    void DynamicComponents::removeEdge(int u, int v) {
        if (G.hasEdge(u, v)) {
            G.removeEdge(u, v);
        }
    }

    bool DynamicComponents::isConnected(int u, int v) const {
        vector<bool> visited(G.numberOfNodes(), false);
        dfs(G, visited, u);
        return visited[v];
    }

    int DynamicComponents::getComponentSize(int v) const {
        vector<bool> visited(G.numberOfNodes(), false);
        dfs(G, visited, v);
        return std::count(visited.begin(), visited.end(), true);
    }

    void dfs(const Graph& G, vector<bool>& visited, int v) {
        visited[v] = true;
        G.forNeighborsOf(v, [&](node u) {
            if (visited[u]) return;
            dfs(G, visited, u);
            });
    }
}

// TODO: find poly-logarithmic version such as:
// TODO: https://github.com/Jonathan-Uy/CSES-Solutions/blob/main/Advanced%20Techniques/Dynamic%20Connectivity.cpp
