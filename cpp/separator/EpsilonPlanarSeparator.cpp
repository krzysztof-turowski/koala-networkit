#include "separator/EpsilonPlanarSeparator.hpp"
#include "networkit/Globals.hpp"
#include "networkit/components/ConnectedComponents.hpp"
#include "separator/GraphUtils.hpp"
#include "separator/PlanarSeparator.hpp"
#include <queue>
#include <vector>

namespace {
double getComponentCost(const std::vector<NetworKit::node> &cc,
                        std::vector<double> &vCost) {
  double c = 0.0;
  for (auto v : cc) {
    c += vCost[v];
  }

  return c;
}

void processSide(
    std::vector<NetworKit::node> partition, std::vector<double> &vertexCost,
    std::unordered_map<int, std::vector<NetworKit::node>> &idToComponentMap,
    double epsilon, std::queue<int> &Q, int &idCnt) {

  double cost = getComponentCost(partition, vertexCost);
  int newId = idCnt++;
  idToComponentMap[newId] = std::move(partition);
  if (cost > epsilon)
    Q.push(newId);
}

} // namespace

namespace Koala {
EpsilonPlanarSeparator::EpsilonPlanarSeparator(
    const NetworKit::Graph &G, double epsilon,
    std::optional<std::map<NetworKit::node, double>> costs)
    : graph(G), epsilon(epsilon) {
  vertexCost.assign(graph.upperNodeIdBound(), 0.0);
  if (costs.has_value()) {
    for (const auto &[v, c] : *costs) {
      vertexCost[v] = c;
    }
  } else {
    NetworKit::count n = graph.numberOfNodes();
    double uniform = n > 0 ? 1.0 / static_cast<double>(n) : 0.0;
    graph.forNodes([&](NetworKit::node v) { vertexCost[v] = uniform; });
  }
}

void EpsilonPlanarSeparator::run() {
  separator.clear();

  int idCnt = 0;
  std::unordered_map<int, std::vector<NetworKit::node>> idToComponentMap;
  std::queue<int> Q;

  NetworKit::ConnectedComponents algoCCs(graph);
  algoCCs.run();

  for (std::vector<NetworKit::node> &cc : algoCCs.getComponents()) {
    double cost = getComponentCost(cc, vertexCost);
    if (cost > epsilon)
      Q.push(idCnt);
    idToComponentMap[idCnt++] = std::move(cc);
  }

  while (!Q.empty()) {
    int curId = Q.front();
    Q.pop();

    auto comp = idToComponentMap[curId];
    auto K = getInducedSubgraph(graph, comp);

    std::vector<double> kCost(comp.size());
    for (size_t i = 0; i < kCost.size(); i++) {
      kCost[i] = vertexCost[comp[i]];
    }

    auto mapFromInducedToOriginal = [&](const std::vector<NetworKit::node> &vec)
        -> std::vector<NetworKit::node> {
      std::vector<NetworKit::node> newVec(vec.size());
      for (size_t i = 0; i < vec.size(); i++) {
        newVec[i] = comp[vec[i]];
      }

      return newVec;
    };

    PlanarSeparator sepAlgo(K, kCost);
    sepAlgo.run();

    for (auto v : sepAlgo.getSeparator())
      separator.push_back(comp[v]);

    processSide(mapFromInducedToOriginal(sepAlgo.getPartitionA()), vertexCost,
                idToComponentMap, epsilon, Q, idCnt);

    processSide(mapFromInducedToOriginal(sepAlgo.getPartitionB()), vertexCost,
                idToComponentMap, epsilon, Q, idCnt);
  }

  hasRun = true;
}
} // namespace Koala
