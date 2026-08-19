#include "separator/MISP.hpp"
#include "networkit/components/ConnectedComponents.hpp"
#include "networkit/graph/GraphTools.hpp"
#include "separator/EpsilonPlanarSeparator.hpp"
#include "separator/GraphUtils.hpp"

namespace Koala {
template <typename T>
MISP<T>::MISP(const NetworKit::Graph &G, double epsilon,
              ComponentSolver componentSolver,
              std::optional<std::map<NetworKit::node, double>> costs)
    : graph(G), epsilon(epsilon), componentSolver(componentSolver),
      costs(costs) {}

template <typename T> void MISP<T>::run() {
  EpsilonPlanarSeparator sep(graph, epsilon, costs);
  sep.run();

  std::unordered_set<NetworKit::node> isInSep(sep.separator.begin(),
                                              sep.separator.end());

  std::unordered_set<NetworKit::node> nonSepNodes;
  graph.forNodes([&](NetworKit::node v) {
    if (isInSep.count(v) == 0) {
      nonSepNodes.insert(v);
    }
  });
  NetworKit::Graph residualGraph =
      NetworKit::GraphTools::subgraphFromNodes(graph, nonSepNodes);

  NetworKit::ConnectedComponents ccAlgo(residualGraph);
  ccAlgo.run();

  for (const auto &componentNodes : ccAlgo.getComponents()) {
    NetworKit::Graph componentGraph =
        getInducedSubgraph(residualGraph, componentNodes);

    for (auto id : componentSolver(componentGraph))
      maximum_independent_set.push_back(componentNodes[id]);
  }
  hasRun = true;
}

template class MISP<NetworKit::node>;
} // namespace Koala
