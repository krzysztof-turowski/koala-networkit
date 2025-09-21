#include <networkit/graph/Graph.hpp>

#pragma once

namespace Koala {

class FlowNetwork {
 public:
  explicit FlowNetwork(const NetworKit::Graph &graph);

  double size() const;

  double lowerCapacity(int u, int v) const;
  double upperCapacity(int u, int v) const;

  void roundFlow();
  void pushValue(int s, int t, double f);

// private:
  const NetworKit::Graph &graph;
  std::vector<std::vector<double>> flow;
  const int N, M;
};

}  // namespace Koala
