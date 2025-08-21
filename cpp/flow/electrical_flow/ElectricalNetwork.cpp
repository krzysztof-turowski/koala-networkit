#include <cassert>
#include <flow/electrical_flow/ElectricalNetwork.hpp>
#include <flow/electrical_flow/LaplaceSolver.hpp>

namespace Koala {
ElectricalNetwork::ElectricalNetwork(const Graph &graph,
                                 const vector<double> &demand)
    : graph(graph), demand(demand) {
  int N = graph.numberOfNodes();
  potentials.assign(N, 0);
  flow.assign(N, vector<double>(N, 0));
};

void ElectricalNetwork::compute(const vector<vector<double>> &resistance) {
  int N = graph.numberOfNodes();
  
  vector<vector<double>> weights(N, vector<double>(N, 0));
  graph.forNodes([&](node u) {
    graph.forNeighborsOf(
        u, [&](node v) { weights[u][v] = 1.0 / resistance[u][v]; });
  });

  VectorXd b(N);
  for (int v = 0; v < N; ++v) {
    b(v) = demand[v];
  }

  auto x = solveLaplace(graph, weights, b);
  for (int v = 0; v < N; ++v) {
    potentials[v] = x(v);
  }
  graph.forNodes([&](node u) {
    graph.forNeighborsOf(u, [&](node v) {
      flow[u][v] = (potentials[u] - potentials[v]) / resistance[u][v];
    });
  });
}

} // namespace Koala
