#include "electric_flow/FlowNetwork.hpp"
#include "networkit/graph/Graph.hpp"

namespace Koala {

FlowNetwork::FlowNetwork(const Graph &graph) : graph(graph) {
  int N = graph.numberOfNodes();

  flow.assign(N, vector<double>(N, 0));
}

double FlowNetwork::upperCapacity(int u, int v) const {
  return graph.weight(u, v) - flow[u][v];
}

double FlowNetwork::lowerCapacity(int u, int v) const {
  return graph.weight(u, v) + flow[u][v];
}

} // namespace Koala
