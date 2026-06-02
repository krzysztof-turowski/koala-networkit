#include <algorithm>
#include <cassert>
#include <cmath>
#include <queue>
#include <vector>

#include "flow/ElectricalFlow.hpp"
#include "flow/electrical_flow/ElectricalNetwork.hpp"
#include "flow/electrical_flow/FlowNetwork.hpp"

namespace Koala {

constexpr double EPSILON = 1e-8;
constexpr double FEASIBILITY_FACTOR = 2.0;
constexpr double ROUTING_COMPLETION_THRESHOLD = 1.0;
constexpr int CONGESTION_NORM = 4;
constexpr double STEP_SIZE_FACTOR = 33.0;
constexpr int VIOLATION_NORM = 2;
constexpr double INTERMEDIATE_VIOLATION_BOUND = 51.0 / 25000.0;
constexpr double FINAL_VIOLATION_BOUND = 1.0 / 100.0;

NetworKit::Graph initialize_graph(
        const NetworKit::Graph &graph, NetworKit::node s, NetworKit::node t) {
    if (!graph.isDirected()) {
        return graph;
    }

    NetworKit::Graph initialized(graph.numberOfNodes(), true, false);
    graph.forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight capacity) {
        if (capacity <= 0 || u == v || u == t || v == s) {
            return;
        }
        initialized.increaseWeight(u, v, capacity);
        initialized.increaseWeight(s, v, capacity);
        initialized.increaseWeight(u, t, capacity);
    });
    initialized.removeMultiEdges();
    return initialized;
}

int get_initial_flow(const NetworKit::Graph &graph, NetworKit::node s, NetworKit::node t) {
    if (!graph.isDirected()) {
        return 0;
    }

    int initial_flow = 0;
    graph.forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight capacity) {
        if (capacity > 0 && u != v && u != t && v != s) {
            initial_flow += capacity;
        }
    });
    return initial_flow;
}

bool push_directed_value(
        const NetworKit::Graph &graph, std::vector<std::vector<double>> &flow,
        NetworKit::node s, NetworKit::node t, double value) {
    NetworKit::count N = graph.numberOfNodes();
    std::vector<std::vector<double>> routed_flow(N, std::vector<double>(N, 0));

    while (value > EPSILON) {
        std::queue<NetworKit::node> queue;
        std::vector<NetworKit::node> parent(N, NetworKit::none);
        std::vector<double> bottleneck(N, 0);
        parent[s] = s;
        bottleneck[s] = value;
        queue.push(s);

        while (!queue.empty() && parent[t] == NetworKit::none) {
            NetworKit::node u = queue.front();
            queue.pop();
            for (NetworKit::node v = 0; v < N; ++v) {
                double capacity = graph.hasEdge(u, v) ? graph.weight(u, v) : 0;
            double residual_capacity = capacity - routed_flow[u][v];
            if (parent[v] == NetworKit::none && residual_capacity > EPSILON) {
                parent[v] = u;
                bottleneck[v] = std::min(bottleneck[u], residual_capacity);
                    queue.push(v);
                }
            }
        }

        if (parent[t] == NetworKit::none) {
            return false;
        }
        double pushed = bottleneck[t];
        for (NetworKit::node v = t; v != s; v = parent[v]) {
            NetworKit::node u = parent[v];
            routed_flow[u][v] += pushed;
            routed_flow[v][u] -= pushed;
        }
        value -= pushed;
    }
    flow.assign(N, std::vector<double>(N, 0));
    for (NetworKit::node u = 0; u < N; ++u) {
        for (NetworKit::node v = 0; v < N; ++v) {
            flow[u][v] = -routed_flow[u][v];
        }
    }
    return true;
}

std::vector<std::vector<double>> get_coupling(const FlowNetwork &f, const std::vector<double> &y) {
    NetworKit::count N = f.graph.numberOfNodes();

    std::vector<std::vector<double>> strength(N, std::vector<double>(N, 0));
    f.graph.forEdges([&](NetworKit::node u, NetworKit::node v) {
        double dy = y[u] - y[v];
        double df = 1.0 / f.upperCapacity(u, v) - 1.0 / f.lowerCapacity(u, v);
        strength[u][v] = dy - df;
        strength[v][u] = df - dy;
    });
    return strength;
}

std::vector<std::vector<double>> get_resistance(const FlowNetwork &f) {
    NetworKit::count N = f.graph.numberOfNodes();

    std::vector<std::vector<double>> resistance(N, std::vector<double>(N, 0));
    f.graph.forEdges([&](NetworKit::node u, NetworKit::node v) {
        resistance[u][v] = resistance[v][u] =
                pow(f.upperCapacity(u, v), -2) + pow(f.lowerCapacity(u, v), -2);
    });
    return resistance;
}

double l_norm(const std::vector<double> &vec, int l) {
    double norm = 0;
    if (l == 0) {
        // infinity
        for (auto x : vec) {
            norm = std::max(norm, x);
        }
        return norm;
    } else if (l > 0) {
        for (auto x : vec) {
            norm += pow(std::abs(x), l);
        }
        return pow(norm, 1.0 / l);
    } else {
        return 0;  // unsupported
    }
}

std::vector<double> get_violation(const FlowNetwork &f, const std::vector<double> &y) {
    std::vector<double> violation(f.graph.numberOfEdges());
    auto coupling = get_coupling(f, y);
    NetworKit::count i = 0;
    f.graph.forEdges([&](NetworKit::node u, NetworKit::node v) {
        violation[i++] = pow(
                coupling[u][v] * std::min(f.upperCapacity(u, v), f.lowerCapacity(u, v)),
                VIOLATION_NORM);
    });
    return violation;
}

ElectricalFlow::ElectricalFlow(
        NetworKit::Graph graph, NetworKit::node s, NetworKit::node t, bool round)
    : MaximumFlow(graph, s, t), original_graph(graph), graph(initialize_graph(graph, s, t)),
      s(s), t(t), U(0),
      initial_flow(get_initial_flow(graph, s, t)), directed(graph.isDirected()), round(round),
            primal(this->graph) {
    this->graph.forNeighborsOf(t, [&](NetworKit::node v) { U += this->graph.weight(v, t); });
}

void ElectricalFlow::run() {
    int L = initial_flow, R = U + 1;
    while (L + 1 < R) {
        target_flow = L + (R - L) / 2;
        if (route_flow()) {
            L = target_flow;
        } else {
            R = target_flow;
        }
    }
    target_flow = L;
    route_flow();

    flow_size = directed ? (L - initial_flow) / 2 : L;
    if (directed) {
        push_directed_value(original_graph, flow, s, t, flow_size);
    } else if (round) {
        primal.flow.assign(graph.numberOfNodes(), std::vector<double>(graph.numberOfNodes(), 0));
        primal.pushValue(s, t, flow_size);
        flow = primal.flow;
    } else {
        flow = primal.flow;
    }
    hasRun = true;
}

bool ElectricalFlow::is_feasible() {
    NetworKit::count M = graph.numberOfEdges();

    double value = 0;
    graph.forNodes([&](NetworKit::node i) { value += demand[i] * dual[i]; });
    return value <= FEASIBILITY_FACTOR * M / (1.0 - progress);
}

bool ElectricalFlow::route_flow() {
    initialize();
    while ((1.0 - progress) * demand[t] > ROUTING_COMPLETION_THRESHOLD) {
        if (!is_feasible()) {
            return false;
        }
        if (!augmentation_step()) {
            return false;
        }
        fixing_step();
    }
    return primal.pushValue(s, t, (1.0 - progress) * demand[t]);
}

void ElectricalFlow::initialize() {
    NetworKit::count N = graph.numberOfNodes();

    demand.assign(N, 0);
    demand[s] = -target_flow, demand[t] = target_flow;
    primal.flow.assign(N, std::vector<double>(N, 0));
    dual.assign(N, 0);
    progress = 0;
}

bool ElectricalFlow::augmentation_step() {
    ElectricalNetwork electrical(graph, demand);
    electrical.compute(get_resistance(primal));

    std::vector<double> congestion(graph.numberOfEdges());
    NetworKit::count i = 0;
    graph.forEdges([&](NetworKit::node u, NetworKit::node v) {
        congestion[i++] = electrical.flow[u][v] / std::min(
                primal.upperCapacity(u, v), primal.lowerCapacity(u, v));
    });
    double congestion_norm = l_norm(congestion, CONGESTION_NORM);
    if (!(congestion_norm > 0) || !std::isfinite(congestion_norm)) {
        return false;
    }
    double step_size = 1.0 / (STEP_SIZE_FACTOR * congestion_norm);

    graph.forEdges([&](NetworKit::node u, NetworKit::node v) {
        primal.flow[u][v] += step_size * electrical.flow[u][v];
        primal.flow[v][u] += step_size * electrical.flow[v][u];
    });
    graph.forNodes([&](NetworKit::node u) { dual[u] += step_size * electrical.potentials[u]; });
    progress += step_size;

    return true;
}

void ElectricalFlow::fixing_step() {
    NetworKit::count N = graph.numberOfNodes();

    auto resistance = get_resistance(primal);
    auto coupling = get_coupling(primal, dual);

    std::vector<std::vector<double>> correction(N, std::vector<double>(N, 0));
    graph.forEdges([&](NetworKit::node u, NetworKit::node v) {
        correction[u][v] = coupling[u][v] / resistance[u][v];
        correction[v][u] = coupling[v][u] / resistance[v][u];
    });

    graph.forEdges([&](NetworKit::node u, NetworKit::node v) {
        primal.flow[u][v] += correction[u][v];
        primal.flow[v][u] += correction[v][u];
    });

    assert(l_norm(get_violation(primal, dual), VIOLATION_NORM) <= INTERMEDIATE_VIOLATION_BOUND);

    std::vector<double> correction_demand(N, 0);
    graph.forEdges([&](NetworKit::node u, NetworKit::node v) {
        correction_demand[u] -= correction[u][v];
        correction_demand[v] -= correction[v][u];
    });

    ElectricalNetwork electrical(graph, correction_demand);
    electrical.compute(get_resistance(primal));

    graph.forEdges([&](NetworKit::node u, NetworKit::node v) {
        primal.flow[u][v] += electrical.flow[u][v];
        primal.flow[v][u] += electrical.flow[v][u];
    });
    graph.forNodes([&](NetworKit::node u) { dual[u] += electrical.potentials[u]; });

    assert(l_norm(get_violation(primal, dual), VIOLATION_NORM) <= FINAL_VIOLATION_BOUND);
}

}  // namespace Koala
