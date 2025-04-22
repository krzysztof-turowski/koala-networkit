#include <networkit/graph/Graph.hpp>
#include <Eigen/Core>

#include <matching/DynamicComponents.hpp>
#include <matching/utils.hpp>

using namespace NetworKit;
using namespace Eigen;
using namespace std;

namespace Koala {
    DynamicComponents::DynamicComponents() {}
    DynamicComponents::DynamicComponents(const Graph& G, const MatrixXd& AG) : G(G), AG(AG) {}

    void DynamicComponents::addVertex(int count) {
        G.addNodes(count);

        MatrixXd AG1 = MatrixXd::Zero(AG.rows() + count, AG.cols() + count);
        AG1.topLeftCorner(AG.rows(), AG.cols()) = AG;
        swap(AG, AG1);
    }

    void DynamicComponents::addVertex() {
        addVertex(1);
    }

    void DynamicComponents::removeVertex(int i) {
        G.removeNode(i);

        int r = AG.rows(), c = AG.cols();
        swap(AG.row(i), AG.row(c-1));
        swap(AG.col(i), AG.col(r-1));
        AG = AG.topLeftCorner(r-1, c-1);
    }

    void DynamicComponents::addEdge(int u, int v) {
        G.addEdge(node(u), node(v));

        AG(u, v) = generateRandom();
        AG(v, u) = -AG(u, v);
    }

    void DynamicComponents::removeEdge(int u, int v) {
        G.removeEdge(node(u), node(v));

        AG(u, v) = 0;
        AG(v, u) = 0;
    }

    int DynamicComponents::size() const {
        return G.numberOfNodes();
    }

    void dfs(const Graph& G, vector<bool>& visited, int v) {
        visited[v] = true;
        G.forNeighborsOf(v, [&](node u) {
            if (visited[u]) return;
            dfs(G, visited, u);
            });
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
}

// TODO: find poly-logarithmic version such as:
// TODO: https://github.com/Jonathan-Uy/CSES-Solutions/blob/main/Advanced%20Techniques/Dynamic%20Connectivity.cpp
