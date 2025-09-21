#include <cassert>
#include <matching/gaussian_matching/BipartiteGaussianMatching.hpp>
#include <matching/gaussian_matching/LazyGaussElimination.hpp>
#include <matching/gaussian_matching/utils.hpp>

#include <NTL/ZZ_p.h>
#include <networkit/graph/DFS.hpp>
#include <networkit/graph/Graph.hpp>

using namespace std;
using namespace NTL;
using namespace NetworKit;

namespace Koala {
pair<vector<int>, vector<int>> getComponents(const Graph &G);

BipartiteGaussianMatching::BipartiteGaussianMatching(const Graph &G1) {
  auto [_G, _oldIdx] = reindexGraph(G1);
  G = _G;
  oldIdx = _oldIdx;

  auto UV = getComponents(G);
  U = UV.first, V = UV.second;

  bpIdx.resize(G.numberOfNodes());
  for (int i = 0; i < U.size(); ++i) {
    bpIdx[U[i]] = i;
  }
  for (int i = 0; i < V.size(); ++i) {
    bpIdx[V[i]] = i;
  }
};

Matching BipartiteGaussianMatching::getMatching() {
  Matching M1;
  for (auto [ui, vi] : M) {
    int u = oldIdx[U[ui]];
    int v = oldIdx[V[vi]];
    M1.insert({u, v});
  }
  return M1;
}

void BipartiteGaussianMatching::run() {
  initZp(ZP_MOD);

  assert(U.size() == V.size());
  int n = max(U.size(), V.size());
  AG = zeroMat(n, n);
  for (auto u : U) {
    for (auto v : G.neighborRange(u)) {
      int ui = bpIdx[u], vi = bpIdx[v];
      auto Xuv = generateRandom();
      AG[ui][vi] = Xuv;
    }
  }

  if (determinant(AG) == 0)
    return;

  MatZp B;
  inv(B, AG);
  auto eliminated = LazyGaussElimination::pivotElimination(
      B, [this](int r, int c) { return AG[r][c] != 0; });

  for (int i = 0; i < eliminated.size(); ++i) {
    M.insert({i, eliminated[i]});
  }
}

pair<vector<int>, vector<int>> getComponents(const Graph &G) {
  int n = G.numberOfNodes();

  vector<int> components[2] = {vector<int>(), vector<int>()};
  vector<int> colors(n, -1);
  bool isBipartite = true;
  for (int r = 0; r < n; ++r) {
    if (colors[r] == -1) {
      colors[r] = 0;
      components[0].push_back(r);

      Traversal::DFSfrom(G, r, [&](node u) {
        G.forNeighborsOf(u, [&](node v) {
          if (colors[v] == -1) {
            colors[v] = (colors[u] + 1) % 2;
            components[colors[v]].push_back(v);
          } else if (colors[u] == colors[v]) {
            isBipartite = false;
          }
        });
      });
    }
  }

  if (!isBipartite) {
    return {};
  }

  return {components[0], components[1]};
}

} // namespace Koala
