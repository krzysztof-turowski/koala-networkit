#include "networkit/graph/Graph.hpp"
#include <vector>

using namespace std;
using namespace NetworKit;

#pragma once
namespace Koala {
class ElectricNetwork {
public:
  ElectricNetwork(const Graph &G, const vector<double> &demand);

  void compute(const vector<vector<double>> &resistance);

  vector<vector<double>> flow;
  vector<double> potentials;

private:
  const Graph &graph;
  const vector<double> &demand;
};
} // namespace Koala
