#include "networkit/Globals.hpp"
#include <flow/electrical_flow/ElectricalFlow.hpp>
#include <flow/electrical_flow/ElectricalNetwork.hpp>
#include <flow/electrical_flow/FlowNetwork.hpp>

using namespace std;
using namespace NetworKit;

namespace Koala {

double currentFlow(const FlowNetwork &f, int t);
vector<vector<double>> getCoupling(const FlowNetwork &f,
                                   const vector<double> &y);
vector<vector<double>> getResistance(const FlowNetwork &f);

void check_augmentation_step(const FlowNetwork &primal,
                             const vector<double> &dual,
                             const vector<double> &violation,
                             const vector<double> &congestion, double stepSize);

double LNorm(const vector<double> &vec, int l) {
  double norm = 0;
  if (l == 0) {
    // infinity
    for (auto x : vec) {
      norm = max(norm, x);
      return norm;
    }
  } else if (l > 0) {
    for (auto x : vec) {
      norm += pow(abs(x), l);
    }
    return pow(norm, 1.0 / l);
  } else {
    return 0;  // unsupported
  }
}

vector<double> getViolation(const FlowNetwork &f, const vector<double> &y) {
  vector<double> violation(f.graph.numberOfEdges());
  auto coupling = getCoupling(f, y);
  int i = 0;
  f.graph.forEdges([&](node u, node v) {
    violation[i++] = pow(
        coupling[u][v] * min(f.lowerCapacity(u, v), f.lowerCapacity(u, v)), 2);
  });
  return violation;
}

ElectricalFlow::ElectricalFlow(const Graph &graph, int s, int t, bool round)
    : graph(graph), s(s), t(t), maximumFlow(0), primal(graph), round(round) {
  U = 0;
  graph.forNeighborsOf(t, [&](node v) { U += graph.weight(v, t); });
}

void ElectricalFlow::run() {
  int L = 0, R = U + 1;
  while (L < R) {
    targetFlow = (L + R) / 2;
    if (routeFlow()) {
      L = targetFlow + 1;
    } else {
      R = targetFlow;
    }
  }
  maximumFlow = L - 1;

  targetFlow = maximumFlow;
  routeFlow();

  if (round) {
    primal.roundFlow();
  }
}

bool ElectricalFlow::isFeasible() {
  int M = graph.numberOfEdges();

  double value = 0;
  graph.forNodes([&](node i) { value += demand[i] * dual[i]; });
  return value <= 2.0 * M / (1.0 - progress);
}

bool ElectricalFlow::routeFlow() {
  initialize();
  while ((1.0 - progress) * demand[t] > 1.0) {
    if (!isFeasible()) {
      return false;
    }
    augmentationStep();
    fixingStep();
  }
  return true;
}

double ElectricalFlow::getFlowSize() const { return maximumFlow; }

void ElectricalFlow::initialize() {
  int N = graph.numberOfNodes();

  demand.assign(N, 0);
  demand[s] = -targetFlow, demand[t] = targetFlow;
  primal.flow.assign(N, vector<double>(N, 0));
  dual.assign(N, 0);
  progress = 0;
}

void ElectricalFlow::augmentationStep() {
  auto violation = getViolation(primal, dual);

  ElectricalNetwork electrical(graph, demand);
  electrical.compute(getResistance(primal));

  vector<double> congestion(graph.numberOfEdges());
  int i = 0;
  graph.forEdges([&](node u, node v) {
    congestion[i++] = electrical.flow[u][v] / min(primal.upperCapacity(u, v),
                                                  primal.lowerCapacity(u, v));
  });
  double stepSize = 1.0 / (33.0 * LNorm(congestion, 4));

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

void ElectricalFlow::fixingStep() {
  int N = graph.numberOfNodes();

  auto resistance = getResistance(primal);
  auto coupling = getCoupling(primal, dual);

  vector<vector<double>> correction(N, vector<double>(N, 0));
  graph.forEdges([&](node u, node v) {
    correction[u][v] = coupling[u][v] / resistance[u][v];
    correction[v][u] = coupling[v][u] / resistance[v][u];
  });

  graph.forEdges([&](node u, node v) {
    primal.flow[u][v] += correction[u][v];
    primal.flow[v][u] += correction[v][u];
  });

  assert(LNorm(getViolation(primal, dual), 2) <= 51.0 / 25000.0);

  vector<double> correctionDemand(N, 0);
  graph.forEdges([&](node u, node v) {
    correctionDemand[u] -= correction[u][v];
    correctionDemand[v] -= correction[v][u];
  });

  ElectricalNetwork electrical(graph, correctionDemand);
  electrical.compute(getResistance(primal));

  graph.forEdges([&](node u, node v) {
    primal.flow[u][v] += electrical.flow[u][v];
    primal.flow[v][u] += electrical.flow[v][u];
  });
  graph.forNodes([&](node u) { dual[u] += electrical.potentials[u]; });

  assert(LNorm(getViolation(primal, dual), 2) <= 1.0 / 100.0);
}

vector<vector<double>> getCoupling(const FlowNetwork &f,
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

vector<vector<double>> getResistance(const FlowNetwork &f) {
  int N = f.graph.numberOfNodes();

  vector<vector<double>> resistance(N, vector<double>(N, 0));
  f.graph.forEdges([&](node u, node v) {
    resistance[u][v] = resistance[v][u] =
        pow(f.upperCapacity(u, v), -2) + pow(f.lowerCapacity(u, v), -2);
  });
  return resistance;
}

double currentFlow(const FlowNetwork &f, int t) {
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
  auto coupling = getCoupling(primal, dual);
  int i = 0;
  primal.graph.forEdges([&](node u, node v) {
    double l = abs(coupling[u][v]) *
               min(primal.upperCapacity(u, v), primal.lowerCapacity(u, v));
    double r =
        4.0 / 3.0 * violation[i] + 7.0 * pow(stepSize * congestion[i], 2);
    assert(l <= r);
    i++;
  });
  assert(LNorm(getViolation(primal, dual), 2) <= 1.0 / 50.0);
}

}  // namespace Koala
