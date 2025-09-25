#include <flow/electrical_flow/ElectricalFlow.hpp>
#include <flow/electrical_flow/ElectricalNetwork.hpp>
#include <flow/electrical_flow/FlowNetwork.hpp>
#include <graph/GraphTools.hpp>
#include <vector>
#include <cmath>
#include <cassert>


using namespace std;
using namespace NetworKit;

namespace Koala {

double current_flow(const FlowNetwork &f, int t);
vector<vector<double>> get_coupling(const FlowNetwork &f,
                                   const vector<double> &y);
vector<vector<double>> get_resistance(const FlowNetwork &f);

void check_augmentation_step(const FlowNetwork &primal,
                             const vector<double> &dual,
                             const vector<double> &violation,
                             const vector<double> &congestion, double stepSize);

double l_norm(const vector<double> &vec, int l) {
  double norm = 0;
  if (l == 0) {
    // infinity
    for (auto x : vec) {
      norm = max(norm, x);
    }
    return norm;
  } else if (l > 0) {
    for (auto x : vec) {
      norm += pow(abs(x), l);
    }
    return pow(norm, 1.0 / l);
  } else {
    return 0;  // unsupported
  }
}

vector<double> get_violation(const FlowNetwork &f, const vector<double> &y) {
  vector<double> violation(f.graph.numberOfEdges());
  auto coupling = get_coupling(f, y);
  int i = 0;
  f.graph.forEdges([&](node u, node v) {
    violation[i++] = pow(
        coupling[u][v] * min(f.lowerCapacity(u, v), f.lowerCapacity(u, v)), 2);
  });
  return violation;
}

ElectricalFlow::ElectricalFlow(Graph graph, int s, int t, bool round)
    : graph(Koala::GraphTools::convertDirectedGraphToUndirected(graph, true)),
      s(s), t(t), U(0), round(round), maximum_flow(0), primal(this->graph) {
  this->graph.forEdges([&](node u, node v) {
    this->graph.setWeight(u, v, this->graph.weight(u, v) / 2.0);
  });
  this->graph.forNeighborsOf(t, [&](node v) { U += this->graph.weight(v, t); });
}

void ElectricalFlow::run() {
  int L = 0, R = U + 1;
  while (L + 1 < R) {
    target_flow = L + (R - L) / 2;
    if (route_flow()) {
      L = target_flow;
    } else {
      R = target_flow;
    }
  }
  target_flow = maximum_flow = L;
  route_flow();

  if (round) {
    primal.roundFlow();
  }
}

bool ElectricalFlow::is_feasible() {
  int M = graph.numberOfEdges();

  double value = 0;
  graph.forNodes([&](node i) { value += demand[i] * dual[i]; });
  return value <= 2.0 * M / (1.0 - progress);
}

bool ElectricalFlow::route_flow() {
  initialize();
  while ((1.0 - progress) * demand[t] > 1.0) {
    if (!is_feasible()) {
      return false;
    }
    augmentation_step();
    fixing_step();
  }
  return true;
}

double ElectricalFlow::getFlowSize() const { return maximum_flow; }

void ElectricalFlow::initialize() {
  int N = graph.numberOfNodes();

  demand.assign(N, 0);
  demand[s] = -target_flow, demand[t] = target_flow;
  primal.flow.assign(N, vector<double>(N, 0));
  dual.assign(N, 0);
  progress = 0;
}

void ElectricalFlow::augmentation_step() {
  auto violation = get_violation(primal, dual);

  ElectricalNetwork electrical(graph, demand);
  electrical.compute(get_resistance(primal));

  vector<double> congestion(graph.numberOfEdges());
  int i = 0;
  graph.forEdges([&](node u, node v) {
    congestion[i++] = electrical.flow[u][v] / min(primal.upperCapacity(u, v),
                                                  primal.lowerCapacity(u, v));
  });
  double stepSize = 1.0 / (33.0 * l_norm(congestion, 4));

  graph.forEdges([&](node u, node v) {
    primal.flow[u][v] += stepSize * electrical.flow[u][v];
    primal.flow[v][u] += stepSize * electrical.flow[v][u];
  });
  graph.forNodes(
      [&](node u) { dual[u] += stepSize * electrical.potentials[u]; });
  progress += stepSize;

  // For debugging purposes
  check_augmentation_step(primal, dual, violation, congestion, stepSize);
}

void ElectricalFlow::fixing_step() {
  int N = graph.numberOfNodes();

  auto resistance = get_resistance(primal);
  auto coupling = get_coupling(primal, dual);

  vector<vector<double>> correction(N, vector<double>(N, 0));
  graph.forEdges([&](node u, node v) {
    correction[u][v] = coupling[u][v] / resistance[u][v];
    correction[v][u] = coupling[v][u] / resistance[v][u];
  });

  graph.forEdges([&](node u, node v) {
    primal.flow[u][v] += correction[u][v];
    primal.flow[v][u] += correction[v][u];
  });

  assert(l_norm(get_violation(primal, dual), 2) <= 51.0 / 25000.0);

  vector<double> correctionDemand(N, 0);
  graph.forEdges([&](node u, node v) {
    correctionDemand[u] -= correction[u][v];
    correctionDemand[v] -= correction[v][u];
  });

  ElectricalNetwork electrical(graph, correctionDemand);
  electrical.compute(get_resistance(primal));

  graph.forEdges([&](node u, node v) {
    primal.flow[u][v] += electrical.flow[u][v];
    primal.flow[v][u] += electrical.flow[v][u];
  });
  graph.forNodes([&](node u) { dual[u] += electrical.potentials[u]; });

  assert(l_norm(get_violation(primal, dual), 2) <= 1.0 / 100.0);
}

vector<vector<double>> get_coupling(const FlowNetwork &f,
                                   const vector<double> &y) {
  int N = f.graph.numberOfNodes();

  vector<vector<double>> strength(N, vector<double>(N, 0));
  f.graph.forEdges([&](node u, node v) {
    double dy = y[u] - y[v];
    double df = 1.0 / f.upperCapacity(u, v) - 1.0 / f.lowerCapacity(u, v);
    strength[u][v] = dy - df;
    strength[v][u] = df - dy;
  });
  return strength;
}

vector<vector<double>> get_resistance(const FlowNetwork &f) {
  int N = f.graph.numberOfNodes();

  vector<vector<double>> resistance(N, vector<double>(N, 0));
  f.graph.forEdges([&](node u, node v) {
    resistance[u][v] = resistance[v][u] =
        pow(f.upperCapacity(u, v), -2) + pow(f.lowerCapacity(u, v), -2);
  });
  return resistance;
}

double current_flow(const FlowNetwork &f, int t) {
  double flow = 0;
  f.graph.forNeighborsOf(t, [&](node v) { flow += f.flow[t][v]; });
  return flow;
}

// Debugging functions
void check_augmentation_step(const FlowNetwork &primal,
                             const vector<double> &dual,
                             const vector<double> &violation,
                             const vector<double> &congestion,
                             double stepSize) {
  auto coupling = get_coupling(primal, dual);
  int i = 0;
  primal.graph.forEdges([&](node u, node v) {
    double l = abs(coupling[u][v]) *
               min(primal.upperCapacity(u, v), primal.lowerCapacity(u, v));
    double r =
        4.0 / 3.0 * violation[i] + 7.0 * pow(stepSize * congestion[i], 2);
    assert(l <= r);
    i++;
  });
  assert(l_norm(get_violation(primal, dual), 2) <= 1.0 / 50.0);
}

}  // namespace Koala
