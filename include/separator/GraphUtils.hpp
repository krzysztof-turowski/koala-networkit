#ifndef KOALA_SEPARATOR_GRAPH_UTILS_HPP
#define KOALA_SEPARATOR_GRAPH_UTILS_HPP

#include "networkit/graph/Graph.hpp"
#include <unordered_map>
#include <vector>

namespace Koala {

inline NetworKit::Graph getInducedSubgraphFromAdj(
    std::vector<std::unordered_set<NetworKit::node>> &adj) {
  NetworKit::Graph induced(adj.size());
  for (size_t i = 0; i < adj.size(); ++i) {
    if (adj[i].size() > 0) {
      for (auto u : adj[i]) {
        if (i < u)
          induced.addEdge(i, u);
      }
    }
  }

  return induced;
}

inline NetworKit::Graph
getInducedSubgraph(const NetworKit::Graph &originalG,
                   const std::vector<NetworKit::node> &nodeIds) {
  std::unordered_map<NetworKit::node, NetworKit::node> idx;
  idx.reserve(nodeIds.size());
  for (std::size_t i = 0; i < nodeIds.size(); ++i)
    idx[nodeIds[i]] = i;

  NetworKit::Graph newG(nodeIds.size());
  for (std::size_t i = 0; i < nodeIds.size(); ++i)
    originalG.forNeighborsOf(nodeIds[i], [&](NetworKit::node w) {
      auto it = idx.find(w);
      if (it != idx.end() && nodeIds[i] < w)
        newG.addEdge(i, it->second);
    });
  return newG;
}

} // namespace Koala

#endif
