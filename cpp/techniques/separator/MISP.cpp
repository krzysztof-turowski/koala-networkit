#include "techniques/separator/MISP.hpp"
#include <algorithm>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "networkit/components/ConnectedComponents.hpp"
#include "networkit/graph/GraphTools.hpp"
#include "techniques/separator/LiptonTarjanPlanarSeparator.hpp"
#include <networkit/graph/AdjListGraph.hpp>

namespace Koala {
template <typename T>
MISP<T>::MISP(const NetworKit::Graph &G, double epsilon, ComponentSolver componentSolver,
              std::optional<std::vector<double>> costs)
    : graph(G), epsilon(epsilon), componentSolver(componentSolver) {
    vertex_cost.assign(G.upperNodeIdBound(), 0.0);
    if (!costs) {
        double average = 1.0 / G.numberOfNodes();
        std::fill(vertex_cost.begin(), vertex_cost.end(), average);
    } else {
        double total_cost = 0.0;
        for (size_t i = 0; i < costs.value().size(); i++) {
            if (G.hasNode(i)) {
                total_cost += costs.value()[i];
                vertex_cost[i] = costs.value()[i];
            }
        }
        assert(total_cost <= 1.0);
    }
}

template <typename T>
double MISP<T>::get_component_cost(const std::vector<NetworKit::node> &connected_component) {
    double c = 0.0;
    for (auto v : connected_component) {
        c += vertex_cost[v];
    }

    return c;
}

template <typename T>
void MISP<T>::processSide(
    const std::vector<NetworKit::node> &subGraph,
    std::unordered_map<NetworKit::node, std::vector<NetworKit::node>> &idToComponentMap,
    std::queue<int> &Q, int &idCnt) {

    double cost = get_component_cost(subGraph);
    int newId = idCnt++;
    idToComponentMap[newId] = std::move(subGraph);
    if (cost > epsilon)
        Q.push(newId);
}

template <typename T>
std::vector<T> MISP<T>::getIndependentSet() {
    assureFinished();
    return maximum_independent_set;
}

template <typename T>
std::vector<NetworKit::node> MISP<T>::find_espilon_planar_separator() {
    std::vector<NetworKit::node> separator;

    int idCnt = 0;
    std::unordered_map<NetworKit::node, std::vector<NetworKit::node>> idToComponentMap;
    std::queue<int> Q;

    NetworKit::ConnectedComponents algoCCs(this->graph);
    algoCCs.run();

    for (std::vector<NetworKit::node> &cc : algoCCs.getComponents()) {
        double cost = get_component_cost(cc);
        if (cost > epsilon)
            Q.push(idCnt);
        idToComponentMap[idCnt++] = std::move(cc);
    }

    while (!Q.empty()) {
        int curId = Q.front();
        Q.pop();

        auto comp = idToComponentMap[curId];
        if (comp.size() <= 5)
            continue;
        std::unordered_set<NetworKit::node> compSet(idToComponentMap[curId].begin(),
                                                    idToComponentMap[curId].end());

        auto induced = NetworKit::GraphTools::subgraphFromNodes(graph, compSet);
        auto kMap = NetworKit::GraphTools::getContinuousNodeIds(induced);
        auto reverseKMap = NetworKit::GraphTools::invertContinuousNodeIds(kMap, induced);
        auto K = NetworKit::GraphTools::getCompactedGraph(induced, kMap);

        std::vector<double> kCost(K.upperNodeIdBound(), 0.0);
        for (auto &[orig, kId] : kMap) {
            kCost[kId] = vertex_cost[orig];
        }
        LiptonTarjanPlanarSeparator sepAlgo(K, kCost);
        sepAlgo.run();

        for (auto v : sepAlgo.getSeparator())
            separator.push_back(reverseKMap[v]);

        auto translatePartition =
            [&](const std::vector<NetworKit::node> &partition) -> std::vector<NetworKit::node> {
            std::vector<NetworKit::node> originalNodes;
            for (auto v : partition) {
                originalNodes.push_back(reverseKMap[v]);
            }
            return originalNodes;
        };

        processSide(translatePartition(sepAlgo.getPartitionA()), idToComponentMap, Q, idCnt);

        processSide(translatePartition(sepAlgo.getPartitionB()), idToComponentMap, Q, idCnt);
    }

    return separator;
}

template <typename T>
void MISP<T>::run() {
    auto separator = find_espilon_planar_separator();

    std::unordered_set<NetworKit::node> isInSep(separator.begin(), separator.end());

    std::unordered_set<NetworKit::node> nonSepNodes;
    graph.forNodes([&](NetworKit::node v) {
        if (isInSep.count(v) == 0) {
            nonSepNodes.insert(v);
        }
    });

    std::vector<NetworKit::node> residualNodesVec(nonSepNodes.begin(), nonSepNodes.end());

    NetworKit::Graph residualGraph = NetworKit::GraphTools::subgraphFromNodes(graph, nonSepNodes);

    NetworKit::ConnectedComponents connectedComponentsAlgorithm(residualGraph);
    connectedComponentsAlgorithm.run();

    for (const auto &componentNodes : connectedComponentsAlgorithm.getComponents()) {
        NetworKit::Graph componentGraph = NetworKit::GraphTools::subgraphFromNodes(
            residualGraph,
            std::unordered_set<NetworKit::node>(componentNodes.begin(), componentNodes.end()));

        auto mapFromOriginalToCompact = NetworKit::GraphTools::getContinuousNodeIds(componentGraph);
        auto mapFromCompactToOriginal = NetworKit::GraphTools::invertContinuousNodeIds(
            mapFromOriginalToCompact, componentGraph);
        auto compactGraph =
            NetworKit::GraphTools::getCompactedGraph(componentGraph, mapFromOriginalToCompact);

        for (auto id : componentSolver(compactGraph)) {
            if constexpr (std::is_same_v<T, NetworKit::node>) {
                maximum_independent_set.push_back(mapFromCompactToOriginal[id]);
            } else if constexpr (std::is_same_v<T, NetworKit::Edge>) {
                maximum_independent_set.push_back(NetworKit::Edge(mapFromCompactToOriginal[id.u],
                                                                  mapFromCompactToOriginal[id.v]));
            }
        }
    }
    hasRun = true;
}
} // namespace Koala
