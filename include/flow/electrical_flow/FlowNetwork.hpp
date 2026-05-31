#pragma once

#include <vector>

#include <networkit/graph/Graph.hpp>

namespace Koala {

class FlowNetwork {
 public:
  explicit FlowNetwork(const NetworKit::Graph &graph);

  double size() const;

  double lowerCapacity(NetworKit::node u, NetworKit::node v) const;
  double upperCapacity(NetworKit::node u, NetworKit::node v) const;

  void roundFlow();
  bool pushValue(NetworKit::node s, NetworKit::node t, double f);

  const NetworKit::Graph &graph;
  std::vector<std::vector<double>> flow;
  const int N, M;
};

}  // namespace Koala
