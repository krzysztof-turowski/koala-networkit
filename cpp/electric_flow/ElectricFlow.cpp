#include <electric_flow/ElectricFlow.hpp>
#include <electric_flow/ElectricNetwork.hpp>
#include <electric_flow/FlowNetwork.hpp>

namespace Koala {

double currentFlow(const FlowNetwork &f, int t);
vector<vector<double>> getCoupling(const FlowNetwork &f,
                                   const vector<double> &y);
vector<vector<double>> getResistance(const FlowNetwork &f);

double LNorm(const vector<double>& vec, int l) {
  double norm = 0;
  if (l == 0) {
    // infinity
    for (auto x: vec) {
      norm = max(norm, x);
      return norm;
    }
  } else if (l>0) {
    for (auto x: vec) {
      norm += pow(abs(x), l);
    }
    return pow(norm, 1.0/l);
  } else {
    return 0; // unsupported
  }
}

vector<double> getViolation(const FlowNetwork &f, const vector<double> &y) {
  vector<double> violation(f.graph.numberOfEdges());
  auto coupling = getCoupling(f, y);
  int i = 0;
  for (auto [u,v]:f.graph.edgeRange()) {
    violation[i++] = pow(coupling[u][v]*min(f.lowerCapacity(u,v), f.lowerCapacity(u,v)),2);
  }
  return violation;
}
double getStepSize(double a, int m);

ElectricFlow::ElectricFlow(const Graph &graph, int s, int t)
    : graph(graph), s(s), t(t), maxFlow(0), primal(graph) {
    }

void ElectricFlow::run() {
  int L = 0, R = primal.U + 1;
  while (L < R) {
    targetFlow = (L + R) / 2;
    if (routeFlow()) {
      L = targetFlow+1;
    } else {
      R = targetFlow;
    }
  }
  maxFlow = L-1;

  targetFlow = maxFlow;
  routeFlow();
  cerr << "remaining:" << targetFlow*(1-progress) << endl;
  // primal.pushValue(s, t, targetFlow*(1-progress));
  primal.roundFlow();
}

bool ElectricFlow::isFeasible() {
  int N = graph.numberOfNodes(), M = graph.numberOfEdges();

  double value = 0;
  for (int i = 0; i < N; ++i) {
    value += demand[i] * dual[i];
  }
  return value <= 2.0 * M / (1.0 - progress);
}

bool ElectricFlow::routeFlow() {
  init();
  while ((1.0 - progress)*demand[t] > 1.0) {
    if (!isFeasible()) {
      return false;
    }
    augmentationStep();
    fixingStep();
  }
  return true;
}

double ElectricFlow::getMaxFlow() const { return maxFlow; }

void ElectricFlow::init() {
  int N = graph.numberOfNodes();

  demand.assign(N, 0);
  demand[s] = -targetFlow, demand[t] = targetFlow;
  primal.flow.assign(N, vector<double>(N, 0));
  dual.assign(N, 0);
  progress = 0;
}

void ElectricFlow::augmentationStep() {
  auto violation = getViolation(primal,dual);

  ElectricNetwork electric(graph, demand);
  electric.compute(getResistance(primal));

  vector<double> congestion(graph.numberOfEdges());
  int i = 0;
  for (auto [u, v]: graph.edgeRange()) {
    congestion[i++] = electric.flow[u][v] / min(primal.upperCapacity(u,v), primal.lowerCapacity(u,v));
  }
  double stepSize = 1/(33.0*LNorm(congestion, 4));

  graph.forNodes([&](node u) {
    graph.forNeighborsOf(
        u, [&](node v) { primal.flow[u][v] += stepSize * electric.flow[u][v]; });
    dual[u] += stepSize * electric.potentials[u];
  });
  progress += stepSize;

  // DEBUG
  {
  auto coupling = getCoupling(primal,dual);
  int i = 0;
  for (auto [u,v]: primal.graph.edgeRange()) {
    bool l = abs(coupling[u][v])*min(primal.upperCapacity(u,v), primal.lowerCapacity(u,v));
    bool r = 4.0/3.0*violation[i] + 7.0*pow(stepSize*congestion[i], 2);
    assert(l <= r);
    i++;
  }
  assert(LNorm(getViolation(primal, dual),2) <= 1.0/50.0);
  }
}

void ElectricFlow::fixingStep() {
  int N = graph.numberOfNodes();

  auto resistance = getResistance(primal);
  auto coupling = getCoupling(primal, dual);

  vector<vector<double>> correction(N, vector<double>(N, 0));
  graph.forNodes([&](node u) {
    graph.forNeighborsOf(u, [&](node v) {
      correction[u][v] = coupling[u][v] / resistance[u][v];
    });
  });

  graph.forNodes([&](node u) {
    graph.forNeighborsOf(
        u, [&](node v) { primal.flow[u][v] += correction[u][v]; });
  });

  assert(LNorm(getViolation(primal, dual),2) <= 51.0/25000.0);

  vector<double> correctionDemand(N, 0);
  for (int u = 0; u < N; ++u) {
    for (int v = 0; v < N; ++v) {
      correctionDemand[u] -= correction[u][v];
    }
  }

  ElectricNetwork electric(graph, correctionDemand);
  electric.compute(getResistance(primal));

  graph.forNodes([&](node u) {
    graph.forNeighborsOf(
        u, [&](node v) { primal.flow[u][v] += electric.flow[u][v]; });
    dual[u] += electric.potentials[u];
  });

  assert(LNorm(getViolation(primal, dual),2) <= 1.0/100.0);
}

vector<vector<double>> getCoupling(const FlowNetwork &f,
                                   const vector<double> &y) {
  int N = f.graph.numberOfNodes();

  vector<vector<double>> strength(N, vector<double>(N, 0));
  f.graph.forNodes([&](node u) {
    f.graph.forNeighborsOf(u, [&](node v) {
      double dy = y[u] - y[v];
      double df = 1.0 / f.upperCapacity(u, v) - 1.0 / f.lowerCapacity(u, v);
      strength[u][v] = dy - df;
    });
  });
  return strength;
}

vector<vector<double>> getResistance(const FlowNetwork &f) {
  int N = f.graph.numberOfNodes();

  vector<vector<double>> resistance(N, vector<double>(N, 0));
  f.graph.forNodes([&](node u) {
    f.graph.forNeighborsOf(u, [&](node v) {
      resistance[u][v] =
          pow(f.upperCapacity(u, v), -2) + pow(f.lowerCapacity(u, v), -2);
    });
  });
  return resistance;
}

double currentFlow(const FlowNetwork &f, int t) {
  double flow = 0;
  for (auto v: f.graph.neighborRange(t)) {
    flow += f.flow[t][v];
  }
  return flow;
}


} // namespace Koala
