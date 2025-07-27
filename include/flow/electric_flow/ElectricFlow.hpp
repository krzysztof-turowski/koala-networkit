#include <flow/electric_flow/FlowNetwork.hpp>
#include <networkit/graph/Graph.hpp>
#include <vector>

using namespace std;
using namespace NetworKit;

#pragma once

namespace Koala {
class ElectricFlow {
public:
  ElectricFlow(const Graph &graph, int s, int t);
  void run();
  double getMaxFlow() const;

// private:
  bool routeFlow();
  void init();
  bool isFeasible();
  void augmentationStep();
  void fixingStep();

  const Graph &graph;
  const int s, t;
  double maxFlow;

  vector<double> demand;
  FlowNetwork primal;
  vector<double> dual;
  double progress;
  double targetFlow;
};

} // namespace Koala
