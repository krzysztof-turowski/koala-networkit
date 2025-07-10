#include "electric_flow/FlowNetwork.hpp"
#include "networkit/graph/Graph.hpp"

namespace Koala {

FlowNetwork::FlowNetwork(const Graph &graph) : graph(graph) {
  int N = graph.numberOfNodes();

  flow.assign(N, vector<double>(N, 0));
}

double FlowNetwork::upperCapacity(int u, int v) const {
  return graph.weight(u, v) + flow[u][v];
}

double FlowNetwork::lowerCapacity(int u, int v) const {
  return graph.weight(u, v) - flow[u][v];
}

// TODO: remove, for debug only
void FlowNetwork::currentDemand() const {
  int N = graph.numberOfNodes();
  vector<double> demand(N, 0);
  graph.forNodes([&](node u){
    graph.forNeighborsOf(u, [&](node v){
      demand[u] += flow[u][v]; 
    });
    cerr << demand[u] << ' ';
  });
  cerr << endl;
}

} // namespace Koala
