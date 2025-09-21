#pragma once

#include <Eigen/Core>
#include <networkit/graph/Graph.hpp>

namespace Koala {
class DynamicComponents {
 public:
  explicit DynamicComponents(const NetworKit::Graph &G);

  void addEdge(int u, int v);
  void removeEdge(int u, int v);
  bool isConnected(int u, int v) const;
  int getComponentSize(int v) const;

  NetworKit::Graph G;
};
}  // namespace Koala
