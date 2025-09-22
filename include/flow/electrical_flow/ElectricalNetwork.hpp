#include <networkit/graph/Graph.hpp>
#include <vector>

#pragma once
namespace Koala {
class ElectricalNetwork {
 public:
  ElectricalNetwork(const NetworKit::Graph &G, const std::vector<double> &demand);

  void compute(const std::vector<std::vector<double>> &resistance);

  std::vector<std::vector<double>> flow;
  std::vector<double> potentials;

 private:
  const NetworKit::Graph &graph;
  const std::vector<double> &demand;
};
}  // namespace Koala
