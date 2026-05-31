#pragma once

#include <vector>

#include <networkit/graph/Graph.hpp>

#include "flow/electrical_flow/FlowNetwork.hpp"

namespace Koala {
class ElectricalFlow {
 public:
  ElectricalFlow(NetworKit::Graph graph, NetworKit::node s, NetworKit::node t,
                 bool round = true);
  void run();
  double getFlowSize() const;

  const NetworKit::Graph& getGraph() const { return originalGraph; }
  const std::vector<std::vector<double>>& getFlow() const { return flow; }

 private:
  bool route_flow();
  void initialize();
  bool is_feasible();
  bool augmentation_step();
  void fixing_step();

  NetworKit::Graph originalGraph;
  NetworKit::Graph graph;
  const NetworKit::node s, t;
  int U;
  int initialFlow;
  bool directed;
  bool round;
  double maximum_flow;

  std::vector<std::vector<double>> flow;
  std::vector<double> demand;
  FlowNetwork primal;
  std::vector<double> dual;
  double progress;
  double target_flow;
};

}  // namespace Koala
