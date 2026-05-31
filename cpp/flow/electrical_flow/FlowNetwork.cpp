#include <algorithm>
#include <utility>
#include <vector>

#include <networkit/graph/Graph.hpp>

#include <flow/electrical_flow/DynamicTree.hpp>
#include <flow/electrical_flow/FlowNetwork.hpp>

double constexpr EPS = 1e-8;

namespace Koala {

FlowNetwork::FlowNetwork(const NetworKit::Graph &graph)
    : graph(graph), N(graph.numberOfNodes()), M(graph.numberOfEdges()) {
  flow.assign(N, std::vector<double>(N, 0));
}

double FlowNetwork::upperCapacity(NetworKit::node u, NetworKit::node v) const {
  return graph.weight(u, v) - flow[u][v];
}

double FlowNetwork::lowerCapacity(NetworKit::node u, NetworKit::node v) const {
  return graph.weight(u, v) + flow[u][v];
}

void cutIntegral(DynamicTree &dt, NetworKit::node u, NetworKit::node v) {
  if (u == v) {
    return;
  }
  auto [mx, my] = dt.pathMin(u, v);
  if (std::abs(dt.weights[mx][my]) <= EPS) {
    dt.cut(mx, my);
    cutIntegral(dt, u, mx);
    cutIntegral(dt, my, v);
  }
}

void FlowNetwork::roundFlow() {
  std::vector<std::vector<double>> weights(N, std::vector<double>(N, 0));
  for (auto [u, v] : graph.edgeRange()) {
    double roundedFlow = std::round(flow[u][v]);
    if (std::abs(flow[u][v] - roundedFlow) <= EPS) {
      flow[u][v] = roundedFlow;
      flow[v][u] = -roundedFlow;
      continue;
    }
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
        std::swap(u, v);
      }
      auto [mx, my] = dt.pathMin(u, v);
      double mf = dt.weights[mx][my] < 0 ? 1.0 + dt.weights[mx][my] : dt.weights[mx][my];
      double vuf = dt.weights[v][u] < 0 ? 1.0 + dt.weights[v][u] : dt.weights[v][u];
      mf = std::min(mf, vuf);

      dt.pathAdd(u, v, mf);
      dt.weights[v][u] -= mf;
      dt.weights[u][v] += mf;

      cutIntegral(dt, u, v);
    }
    if (std::abs(dt.weights[u][v]) > EPS) {
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

bool FlowNetwork::pushValue(NetworKit::node s, NetworKit::node t, double f) {
  while (f > EPS) {
    std::vector<NetworKit::node> st;
    std::vector<std::pair<NetworKit::node, double>> parent(
        graph.numberOfNodes(), {NetworKit::none, 0.0});

    parent[t] = {t, f};
    st.push_back(t);

    while (!st.empty()) {
      NetworKit::node v = st.back();
      st.pop_back();
      if (v == s) {
        break;
      }

      graph.forNeighborsOf(v, [&](NetworKit::node u) {
        double residualCapacity = graph.weight(u, v) - flow[v][u];
        if (parent[u].first == NetworKit::none && residualCapacity > EPS) {
          st.push_back(u);
          parent[u] = {v, std::min(parent[v].second, residualCapacity)};
        }
      });
    }

    if (parent[s].first == NetworKit::none) {
      return false;
    }
    NetworKit::node v = s;
    double f1 = parent[s].second;
    while (v != t) {
      flow[parent[v].first][v] += f1;
      flow[v][parent[v].first] -= f1;
      v = parent[v].first;
    }
    f -= f1;
  }
  return true;
}
}  // namespace Koala
