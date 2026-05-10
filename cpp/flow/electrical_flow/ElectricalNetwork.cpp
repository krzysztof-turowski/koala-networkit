#include <cassert>
#include <vector>

#include <networkit/graph/Graph.hpp>

#include <flow/electrical_flow/ElectricalNetwork.hpp>
#include <flow/electrical_flow/LaplaceSolver.hpp>

namespace Koala {

ElectricalNetwork::ElectricalNetwork(
    const NetworKit::Graph &graph, const std::vector<double> &demand)
        : graph(graph), demand(demand) {
  int N = graph.numberOfNodes();
  potentials.assign(N, 0);
  flow.assign(N, std::vector<double>(N, 0));
}

void ElectricalNetwork::compute(const std::vector<std::vector<double>> &resistance) {
  int N = graph.numberOfNodes();

  std::vector<std::vector<double>> weights(N, std::vector<double>(N, 0));
  graph.forEdges([&](NetworKit::node u, NetworKit::node v) {
    weights[u][v] = 1.0 / resistance[u][v];
    weights[v][u] = 1.0 / resistance[v][u];
  });

  Eigen::VectorXd b(N);
  for (int v = 0; v < N; ++v) {
    b(v) = demand[v];
  }

  auto x = solveLaplace(graph, weights, b);
  for (int v = 0; v < N; ++v) {
    potentials[v] = x(v);
  }
  graph.forEdges([&](NetworKit::node u, NetworKit::node v) {
    flow[u][v] = (potentials[u] - potentials[v]) / resistance[u][v];
    flow[v][u] = -flow[u][v];
  });
}

}  // namespace Koala
