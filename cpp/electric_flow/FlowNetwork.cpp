#include "networkit/graph/Graph.hpp"
#include "electric_flow/FlowNetwork.hpp"
#include "electric_flow/DynamicTree.hpp"

double constexpr EPS = 1e-8;

namespace Koala {

inline int getU(const Graph& graph) {
  int U = 0;
  for (auto [u,v]: graph.edgeRange()) {
    U += abs(graph.weight(u,v));
  }
  return U;
}

FlowNetwork::FlowNetwork(const Graph &graph) : graph(graph), N(graph.numberOfNodes()), M(graph.numberOfEdges()), U(getU(graph)) {
  flow.assign(N, vector<double>(N, 0));
}

double FlowNetwork::upperCapacity(int u, int v) const {
  return graph.weight(u, v) - flow[u][v];
}

double FlowNetwork::lowerCapacity(int u, int v) const {
  return graph.weight(u, v) + flow[u][v];
}

void FlowNetwork::roundFlow() {
  vector<vector<double>> weights(N, vector<double>(N, 0));
  for (auto [u, v]: graph.edgeRange()) {
    double intFlow;
    weights[u][v] = modf(flow[u][v], &intFlow);
    weights[v][u] = -weights[u][v];

    flow[u][v] = intFlow;
    flow[v][u] = -intFlow;
  }

  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      cout << weights[i][j] << ' ';
    }
    cout << endl;
  }
  cout << endl;

  DynamicTree dt(N, weights);

  for (auto [u,v]: graph.edgeRange()) {
    if (dt.findRoot(u) == dt.findRoot(v)) {
      double s = dt.pathSum(u,v) + dt.weights[v][u];
      if (s < 0) {
        swap(u,v);
      }
      auto [mx,my] = dt.pathMin(u,v);
      if (dt.weights[v][u] >= 0 && dt.weights[v][u] < dt.weights[mx][my]) {
        mx=v, my=u;
      }
      double mf = dt.weights[mx][my];
      dt.pathAdd(u, v, mf);
      dt.weights[v][u] -= mf;
      dt.weights[u][v] += mf;

      for (auto [x,y]: dt.graph.edgeRange()) {
        if (abs(dt.weights[x][y]) <= EPS) {
          dt.cut(x, y);
        }
      }
    }
    if (abs(dt.weights[u][v]) > EPS) {
      dt.link(u,v);
    }
  }

  // for (auto [u,v]: dt.graph.edgeRange()) {
    // dt.weights[v][v] = round(dt.weights[v][v]);
    // dt.weights[v][u] = round(dt.weights[v][v]);
    // cout << u << ' ' << v << ' ' << dt.weights[u][v] << endl;
  // }

  for (int  i = 0; i < N; ++i) {
  for (int  j = 0; j < N; ++j)
    cout << dt.weights[i][j] << ' ';
  cout << endl;
  }


  // for (auto [u,v]: graph.edgeRange()) {
  //     flow[u][v] += dt.weights[u][v];
  //     flow[v][u] -= dt.weights[u][v];
  // }
}

void FlowNetwork::pushValue(int s, int t, double f) {\
  while (f > EPS) {
    vector<int> st;
    vector<pair<int, double>> parent(graph.numberOfNodes(), {-1, 0.0});

    parent[t] = {t, f};
    st.push_back(t);

    while (!st.empty()) {
      int v = st.back();
      st.pop_back();

      for (auto u: graph.neighborRange(v)) {
        if (parent[u].first == -1) {
          st.push_back(u);
          parent[u] = {v, min(parent[v].second, graph.weight(u,v)-flow[u][v])};
        }
      }
    }

    int v = s;
    double f1 = parent[s].second;
    while (v != t) {
      flow[parent[v].first][v] += f1;
      flow[v][parent[v].first] -= f1;
      v = parent[v].first;
    }
    f-=f1;
  }
}
} // namespace Koala
