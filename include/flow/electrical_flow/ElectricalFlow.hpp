#include <flow/electrical_flow/FlowNetwork.hpp>
#include <networkit/graph/Graph.hpp>
#include <vector>

#pragma once

namespace Koala {
class ElectricalFlow {
public:
  ElectricalFlow(const NetworKit::Graph &graph, int s, int t, bool round=true);
  void run();
  double getFlowSize() const;

// private:
  bool routeFlow();
  void initialize();
  bool isFeasible();
  void augmentationStep();
  void fixingStep();

  const NetworKit::Graph &graph;
  const int s, t;
  int U;
  bool round;
  double maximumFlow;

  std::vector<double> demand;
  FlowNetwork primal;
  std::vector<double> dual;
  double progress;
  double targetFlow;
};

} // namespace Koala
