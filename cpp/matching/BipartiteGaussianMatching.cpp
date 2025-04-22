#include <matching/GaussElimination.hpp>
#include <matching/BipartiteGaussianMatching.hpp>
#include <matching/utils.hpp>

#include <networkit/graph/Graph.hpp>
#include <Eigen/Core>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;
using namespace NetworKit;

namespace Koala {
    bool dfs(const Graph& G, int v, vector<int>& colors);
    pair<vector<int>, vector<int>> getComponents(const Graph& G);

    BipartiteGaussianMatching::BipartiteGaussianMatching(const Graph& G) : G(G) {
        auto UV = getComponents(G);

        U = UV.first, V = UV.second;
        newIdx.resize(G.numberOfNodes());
        for (int i = 0; i < U.size(); ++i) {
            newIdx[U[i]] = i;
        }
        for (int i = 0; i < V.size(); ++i) {
            newIdx[V[i]] = i;
        }
    };

    Matching BipartiteGaussianMatching::getMatching() {
        Matching M1;
        for (auto [ui, vi]: M) {
            M1.insert({U[ui], V[vi]});
        }
        return M1;
    }

    void BipartiteGaussianMatching::run() {
        int n = max(U.size(), V.size());
        AG = ArrayXXd::Zero(n, n);
        for (auto& u : U) {
            G.forEdgesOf(u, [&](int v) {
                int ui = newIdx[u], vi = newIdx[v];
                double Xuv = generateRandom();
                AG(ui, vi) = Xuv;
                });
        }
        
        if (eq(AG.determinant(), 0))
            return;

        MatrixXd B = AG.inverse();
        auto eliminated = pivotElimination(B, [this](int r, int c) {return !eq(AG(r, c), 0);});

        for (int i = 0; i < eliminated.size(); ++i) {
            M.insert({ i, eliminated[i] });
        }
    }

    pair<vector<int>, vector<int>> getComponents(const Graph& G) {
        int n = G.numberOfNodes();

        vector<int> colors(n, -1);
        for (int v = 0; v < n; ++v) {
            if (colors[v] == -1) {
                colors[v] = 0;
                assert(!dfs(G, v, colors));
            }
        }

        vector<int> res[2] = { vector<int>(), vector<int>() };
        for (int v = 0; v < n; ++v) {
            res[colors[v]].push_back(v);
        }
        return { res[0], res[1] };
    }

    bool dfs(const Graph& G, int v, vector<int>& colors) {
        bool broken = false;
        G.forEdgesOf(v, [&](int u) {
            if (broken) return;
            if (colors[u] == -1) {
                colors[u] = (colors[v] + 1) % 2;
                broken = dfs(G, u, colors);
            } else if (colors[v] == colors[u]) {
                broken = true;
            }
            });
        return broken;
    };
}
