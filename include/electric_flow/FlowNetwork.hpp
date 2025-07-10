#include "networkit/graph/Graph.hpp"

using namespace std;
using namespace NetworKit;

#pragma once

namespace Koala {

class FlowNetwork {
public:
  FlowNetwork(const Graph &graph);

  double size() const;

  double lowerCapacity(int u, int v) const;
  double upperCapacity(int u, int v) const;
  void currentDemand() const;

// private:
  const Graph &graph;
  vector<vector<double>> flow;
};

} // namespace Koala
