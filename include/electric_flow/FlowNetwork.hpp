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

  void roundFlow();
  void pushValue(int s, int t, double f);

// private:
  const Graph &graph;
  vector<vector<double>> flow;
  const int N, M, U;
};

} // namespace Koala
