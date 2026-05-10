#pragma once

#include <networkit/graph/Graph.hpp>

namespace Koala {
class DynamicComponents {
 public:
  explicit DynamicComponents(const NetworKit::Graph &G);

  void addEdge(NetworKit::node u, NetworKit::node v);
  void removeEdge(NetworKit::node u, NetworKit::node v);
  bool isConnected(NetworKit::node u, NetworKit::node v) const;
  int getComponentSize(NetworKit::node v) const;

  NetworKit::Graph G;
};
}  // namespace Koala
