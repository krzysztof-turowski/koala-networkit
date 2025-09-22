#include "networkit/graph/Graph.hpp"
#include <cassert>
#include <flow/electrical_flow/ElectricalNetwork.hpp>
#include <flow/electrical_flow/LaplaceSolver.hpp>

using namespace std;
using namespace NetworKit;
using namespace Eigen;

namespace Koala {

ElectricalNetwork::ElectricalNetwork(const Graph &graph,
                                     const vector<double> &demand)
    : graph(graph), demand(demand) {
  int N = graph.numberOfNodes();
  potentials.assign(N, 0);
  flow.assign(N, vector<double>(N, 0));
}

void ElectricalNetwork::compute(const vector<vector<double>> &resistance) {
  int N = graph.numberOfNodes();

  vector<vector<double>> weights(N, vector<double>(N, 0));
  graph.forEdges([&](node u, node v) {
    weights[u][v] = 1.0 / resistance[u][v];
    weights[v][u] = 1.0 / resistance[v][u];
  });

  VectorXd b(N);
  for (int v = 0; v < N; ++v) {
    b(v) = demand[v];
  }

  auto x = solveLaplace(graph, weights, b);
  for (int v = 0; v < N; ++v) {
    potentials[v] = x(v);
  }
  graph.forEdges([&](node u, node v) {
    flow[u][v] = (potentials[u] - potentials[v]) / resistance[u][v];
    flow[v][u] = -flow[u][v];
  });
}

}  // namespace Koala
