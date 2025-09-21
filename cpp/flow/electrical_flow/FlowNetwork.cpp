#include "flow/electrical_flow/FlowNetwork.hpp"
#include "flow/electrical_flow/DynamicTree.hpp"
#include "networkit/graph/Graph.hpp"

double constexpr EPS = 1e-8;

using namespace std;
using namespace NetworKit;

namespace Koala {

FlowNetwork::FlowNetwork(const Graph &graph)
    : graph(graph), N(graph.numberOfNodes()), M(graph.numberOfEdges()) {
  flow.assign(N, vector<double>(N, 0));
}

double FlowNetwork::upperCapacity(int u, int v) const {
  return graph.weight(u, v) - flow[u][v];
}

double FlowNetwork::lowerCapacity(int u, int v) const {
  return graph.weight(u, v) + flow[u][v];
}

void cutIntegral(DynamicTree &dt, int u, int v) {
  if (u == v) {
    return;
  }
  auto [mx, my] = dt.pathMin(u, v);
  if (abs(dt.weights[mx][my]) <= EPS) {
    dt.cut(mx, my);
    cutIntegral(dt, u, mx);
    cutIntegral(dt, my, v);
  }
}

void FlowNetwork::roundFlow() {
  vector<vector<double>> weights(N, vector<double>(N, 0));
  for (auto [u, v] : graph.edgeRange()) {
    double intFlow;
    weights[u][v] = modf(flow[u][v], &intFlow);
    weights[v][u] = -weights[u][v];

    flow[u][v] = intFlow;
    flow[v][u] = -intFlow;
  }

  DynamicTree dt(N, weights);

  for (auto [u, v] : graph.edgeRange()) {
    if (dt.findRoot(u) == dt.findRoot(v)) {
      double s = dt.pathSum(u, v) + dt.weights[v][u];
      if (s < 0) {
        swap(u, v);
      }
      auto [mx, my] = dt.pathMin(u, v);
      double mf = dt.weights[mx][my] < 0 ? 1.0 + dt.weights[mx][my]
                                         : dt.weights[mx][my];
      double vuf =
          dt.weights[v][u] < 0 ? 1.0 + dt.weights[v][u] : dt.weights[v][u];
      mf = min(mf, vuf);

      dt.pathAdd(u, v, mf);
      dt.weights[v][u] -= mf;
      dt.weights[u][v] += mf;

      cutIntegral(dt, u, v);
    }
    if (abs(dt.weights[u][v]) > EPS) {
      dt.link(u, v);
    }
  }

  for (auto [u, v] : graph.edgeRange()) {
    double f = dt.weights[u][v] >= 0 ? ceil(dt.weights[u][v])
                                     : floor(dt.weights[u][v]);
    flow[u][v] += f;
    flow[v][u] -= f;
  }
}

void FlowNetwork::pushValue(int s, int t, double f) {
  while (f > EPS) {
    vector<int> st;
    vector<pair<int, double>> parent(graph.numberOfNodes(), {-1, 0.0});

    parent[t] = {t, f};
    st.push_back(t);

    while (!st.empty()) {
      int v = st.back();
      st.pop_back();

      graph.forNeighborsOf(v, [&](node u) {
        if (parent[u].first == -1) {
          st.push_back(u);
          parent[u] = {v,
                       min(parent[v].second, graph.weight(u, v) - flow[u][v])};
        }
      });
    }

    int v = s;
    double f1 = parent[s].second;
    while (v != t) {
      flow[parent[v].first][v] += f1;
      flow[v][parent[v].first] -= f1;
      v = parent[v].first;
    }
    f -= f1;
  }
}
}  // namespace Koala
