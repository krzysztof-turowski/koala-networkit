#include <networkit/graph/Graph.hpp>

#include "matching/gaussian_matching/DynamicComponents.hpp"
#include "matching/gaussian_matching/utils.hpp"

namespace Koala {
  void dfs(const NetworKit::Graph &G, std::vector<bool> &visited, int v) {
    visited[v] = true;
    G.forNeighborsOf(v, [&](NetworKit::node u) {
      if (visited[u])
        return;
      dfs(G, visited, u);
    });
  }

  DynamicComponents::DynamicComponents(const NetworKit::Graph &G) : G(G) {}

  void DynamicComponents::addEdge(NetworKit::node u, NetworKit::node v) { G.addEdge(u, v); }

  void DynamicComponents::removeEdge(NetworKit::node u, NetworKit::node v) {
    if (G.hasEdge(u, v)) {
      G.removeEdge(u, v);
    }
}

bool DynamicComponents::isConnected(NetworKit::node u, NetworKit::node v) const {
  std::vector<bool> visited(G.numberOfNodes(), false);
  dfs(G, visited, u);
  return visited[v];
}

int DynamicComponents::getComponentSize(NetworKit::node v) const {
  std::vector<bool> visited(G.numberOfNodes(), false);
  dfs(G, visited, v);
  return std::count(visited.begin(), visited.end(), true);
}
}  // namespace Koala

