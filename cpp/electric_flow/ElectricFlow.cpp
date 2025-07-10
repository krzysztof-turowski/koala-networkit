#include <electric_flow/ElectricFlow.hpp>
#include <electric_flow/ElectricNetwork.hpp>
#include <electric_flow/FlowNetwork.hpp>

double constexpr STEP = 0.1;

namespace Koala {
vector<vector<double>> getCoupling(const FlowNetwork &f,
                                   const vector<double> &y);
vector<vector<double>> getResistance(const FlowNetwork &f);

ElectricFlow::ElectricFlow(const Graph &graph, int s, int t)
    : graph(graph), s(s), t(t), maxFlow(0), primal(graph) {}

void ElectricFlow::run() {
  cerr << "RUN" << endl;
  /*
    double L = 0, R = 1000;
    while (L < R) {
      double M = (L + R) / 2;
      if (routeFlow(M)) {
        L = M;
      } else {
        R = M;
      }
    }
    maxFlow = L;*/
  targetFlow = 10;
  maxFlow = routeFlow() ? 10 : 0;
}

bool ElectricFlow::isFeasable() {
  int N = graph.numberOfNodes(), M = graph.numberOfEdges();

  double value = 0;
  for (int i = 0; i < N; ++i) {
    value += demand[i] * dual[i];
  }
  return value <= 2.0 * M / (1.0 - progress);
}

bool ElectricFlow::routeFlow() {
  cerr << "ROUTE_FLOW: " << targetFlow << endl;

  init();
  int it = 1.0 / STEP;
  for (int i = 0; i < it; ++i) {
    augmentationStep();
    fixingStep();

    if (!isFeasable()) {
      return false;
    }
    cerr << "PROGRES: " << progress << endl;
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
  cerr << "AUGMENTATION" << endl;

  ElectricNetwork electric(graph, demand);
  electric.compute(getResistance(primal));

  graph.forNodes([&](node u) {
    graph.forNeighborsOf(
        u, [&](node v) { primal.flow[u][v] += STEP * electric.flow[u][v]; });
    dual[u] += STEP * electric.potentials[u];
  });
  progress += STEP;
}

void ElectricFlow::fixingStep() {
  cerr << "FIXING" << endl;

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
}

vector<vector<double>> getCoupling(const FlowNetwork &f,
                                   const vector<double> &y) {
  int N = f.graph.numberOfNodes();

  vector<vector<double>> strength(N, vector<double>(N, 0));
  f.graph.forNodes([&](node u) {
    f.graph.forNeighborsOf(u, [&](node v) {
      double dy = y[v] - y[u];
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

} // namespace Koala
