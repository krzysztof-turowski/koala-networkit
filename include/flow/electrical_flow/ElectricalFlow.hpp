#pragma once

#include <vector>

#include <networkit/graph/Graph.hpp>

#include "flow/electrical_flow/FlowNetwork.hpp"

class GenTest;

namespace Koala {
class ElectricalFlow {
 public:
  ElectricalFlow(NetworKit::Graph graph, int s, int t,
                 bool round = true);
  void run();
  double getFlowSize() const;

  const NetworKit::Graph& getGraph() const { return graph; }
  const std::vector<std::vector<double>>& getFlow() const { return primal.flow; }

 private:
  bool route_flow();
  void initialize();
  bool is_feasible();
  void augmentation_step();
  void fixing_step();

  NetworKit::Graph graph;
  const int s, t;
  int U;
  bool round;
  double maximum_flow;

  std::vector<double> demand;
  FlowNetwork primal;
  std::vector<double> dual;
  double progress;
  double target_flow;
};

}  // namespace Koala
