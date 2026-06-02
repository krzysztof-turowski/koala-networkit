#include <algorithm>
#include <utility>
#include <vector>

#include <networkit/graph/Graph.hpp>

#include "flow/electrical_flow/FlowNetwork.hpp"
#include "structures/dynamic_tree/NaiveDynamicTree.hpp"

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

void cut_integral(
        NaiveDynamicTree<double> &dt, NetworKit::node u, NetworKit::node v) {
    if (u == v) {
        return;
    }
    auto [mx, my] = dt.pathMin(u, v);
    if (std::abs(dt.getWeight(mx, my)) <= EPS) {
        dt.cut(mx, my);
        cut_integral(dt, u, mx);
        cut_integral(dt, my, v);
    }
}

void FlowNetwork::roundFlow() {
    std::vector<std::vector<double>> weights(N, std::vector<double>(N, 0));
    for (auto [u, v] : graph.edgeRange()) {
        double rounded_flow = std::round(flow[u][v]);
        if (std::abs(flow[u][v] - rounded_flow) <= EPS) {
            flow[u][v] = rounded_flow;
            flow[v][u] = -rounded_flow;
            continue;
        }
        double int_flow;
        weights[u][v] = modf(flow[u][v], &int_flow);
        weights[v][u] = -weights[u][v];

        flow[u][v] = int_flow;
        flow[v][u] = -int_flow;
    }

    NaiveDynamicTree<double> dt(N, weights);

    for (auto [u, v] : graph.edgeRange()) {
        if (dt.findRoot(u) == dt.findRoot(v)) {
            double s = dt.pathSum(u, v) + dt.getWeight(v, u);
            if (s < 0) {
                std::swap(u, v);
            }
            auto [mx, my] = dt.pathMin(u, v);
            double mf = dt.getWeight(mx, my) < 0 ? 1.0 + dt.getWeight(mx, my)
                                                 : dt.getWeight(mx, my);
            double vuf = dt.getWeight(v, u) < 0 ? 1.0 + dt.getWeight(v, u)
                                               : dt.getWeight(v, u);
            mf = std::min(mf, vuf);

            dt.pathAdd(u, v, mf);
            dt.addWeight(v, u, -mf);
            dt.addWeight(u, v, mf);

            cut_integral(dt, u, v);
        }
        if (std::abs(dt.getWeight(u, v)) > EPS) {
            dt.link(u, v, dt.getWeight(u, v));
        }
    }

    for (auto [u, v] : graph.edgeRange()) {
        double f = dt.getWeight(u, v) >= 0 ? ceil(dt.getWeight(u, v))
                                          : floor(dt.getWeight(u, v));
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
                double residual_capacity = graph.weight(u, v) - flow[v][u];
                if (parent[u].first == NetworKit::none && residual_capacity > EPS) {
                    st.push_back(u);
                    parent[u] = {v, std::min(parent[v].second, residual_capacity)};
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
