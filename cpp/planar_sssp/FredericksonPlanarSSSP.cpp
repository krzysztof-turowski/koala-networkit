#include "planar_sssp/FredericksonPlanarSSSP.hpp"

#include <algorithm>
#include <cmath>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/distance/Dijkstra.hpp>
#include <networkit/distance/MultiTargetDijkstra.hpp>
#include <networkit/graph/BFS.hpp>
#include <networkit/graph/DFS.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <utility>

#include "planar_sssp/PlanarUtils.hpp"
#include "planar_sssp/SuitableRDivision.hpp"
#include "planar_sssp/TopologyBasedHeap.hpp"

using pairDistance_t =
    std::unordered_map<std::pair<int, int>, int, boost::hash<std::pair<int, int>>>;
using nodeSubsets_t = std::vector<std::vector<NetworKit::node>>;

namespace Koala {
pairDistance_t getDistancebetweenBoundaryNodesLevel2(NetworKit::Graph& graph,
    nodeSubsets_t& division, const std::unordered_set<NetworKit::node>& extraBoundaryNodes) {
    std::vector<pairDistance_t> resultPerRegion(division.size());
    pairDistance_t result;

    std::vector<int> regionCount(graph.numberOfNodes(), 0);
    std::vector<int> isBoundary(graph.numberOfNodes(), 0);

    for (int i = 0; i < division.size(); i++) {
        for (auto node : division[i]) {
            regionCount[node]++;
        }
    }
    for (int i = 0; i < graph.numberOfNodes(); i++) {
        if (regionCount[i] > 1) {
            isBoundary[i] = true;
        }
    }

    for (int r = 0; r < division.size(); r++) {
        auto& region = division[r];
        auto subGraph = NetworKit::GraphTools::subgraphFromNodes(
            graph, std::unordered_set<NetworKit::node>{region.begin(), region.end()});
        std::vector<NetworKit::node> boundary;

        for (auto node : region) {
            if (isBoundary[node] || (extraBoundaryNodes.find(node) != extraBoundaryNodes.end())) {
                boundary.push_back(node);
            }
        }

        for (int i = 0; i < boundary.size(); i++) {
            NetworKit::MultiTargetDijkstra dijkstra(
                subGraph, boundary[i], boundary.begin(), boundary.end());
            dijkstra.run();
            auto map = dijkstra.getTargetIndexMap();
            auto distances = dijkstra.getDistances();
            for (int j = 0; j < boundary.size(); j++) {
                if (distances[map[boundary[j]]] < 0 ||
                    distances[map[boundary[j]]] > std::numeric_limits<int>::max())
                    continue;
                resultPerRegion[r][std::make_pair(boundary[i], boundary[j])] =
                    distances[map[boundary[j]]];
                resultPerRegion[r][std::make_pair(boundary[j], boundary[i])] =
                    distances[map[boundary[j]]];
            }
        }
    }

    for (auto pairs : resultPerRegion) {
        for (auto [pair, distance] : pairs) {
            if (result.contains(pair)) {
                result[pair] = std::min(result[pair], distance);
            } else {
                result[pair] = distance;
            }
        }
    }

    return result;
}

int mopUp(std::vector<int>& shortestDistance, nodeSubsets_t& regions, NetworKit::Graph& G,
    int targetRegion, int t) {
    int result = std::numeric_limits<int>::max();

    std::unordered_set<NetworKit::node> subGraphNodes(
        regions[targetRegion].begin(), regions[targetRegion].end());
    NetworKit::Graph subGraph = NetworKit::GraphTools::subgraphFromNodes(G, subGraphNodes);

    NetworKit::MultiTargetDijkstra dijkstra(
        subGraph, t, regions[targetRegion].begin(), regions[targetRegion].end());
    dijkstra.run();
    auto map = dijkstra.getTargetIndexMap();
    auto distances = dijkstra.getDistances();

    for (auto node : regions[targetRegion]) {
        if (shortestDistance[node] == -1) continue;
        result = std::min(result, shortestDistance[node] + static_cast<int>(distances[map[node]]));
    }
    return result;
}

std::vector<int> mainThrust(NetworKit::Graph& graph, nodeSubsets_t& regions,
    pairDistance_t& distances, int s, const std::unordered_set<NetworKit::node>& extraBoudryNodes) {
    std::vector<int> shortestDistance(graph.numberOfNodes(), -1);

    TopologyHeap heap(graph, regions, distances, s, extraBoudryNodes);

    while (!heap.empty()) {
        auto [distance, node] = heap.top();
        if (node == -1) break;
        shortestDistance[node] = distance;
        heap.closeNode(node);
    }

    return shortestDistance;
}

pairDistance_t getDistanceInRegionLevel1(
    NetworKit::Graph& graph, std::vector<NetworKit::node>& level1BoundaryNodes, int r2, int c) {
    pairDistance_t result;

    NetworKit::ConnectedComponents cc(graph);
    cc.run();
    nodeSubsets_t components = cc.getComponents();

    for (auto& component : components) {
        auto connectedSubgraph =
            NetworKit::GraphTools::subgraphFromNodes(graph, {component.begin(), component.end()});
        std::vector<NetworKit::node> nodes;
        std::unordered_map<NetworKit::node, NetworKit::node> mapOfNodes;

        for (auto node : connectedSubgraph.nodeRange()) nodes.push_back(node);
        for (int i = 0; i < nodes.size(); i++) mapOfNodes[nodes[i]] = i;

        NetworKit::Graph subGraphFrom0(nodes.size(), true);
        for (auto e : connectedSubgraph.edgeWeightRange()) {
            subGraphFrom0.setWeight(mapOfNodes[e.u], mapOfNodes[e.v], e.weight);
        }

        nodeSubsets_t level2Division;
        if (subGraphFrom0.numberOfNodes() <= r2) {
            level2Division.push_back(std::vector<NetworKit::node>{});
            for (auto node : subGraphFrom0.nodeRange()) {
                level2Division[0].push_back(node);
            }
        } else {
            level2Division = findSuitableRDivision(subGraphFrom0, r2, c);
        }

        std::unordered_set<NetworKit::node> boundaryNodesOfComponent;
        for (auto bn : level1BoundaryNodes) {
            if (mapOfNodes.find(bn) != mapOfNodes.end()) {
                boundaryNodesOfComponent.insert(mapOfNodes[bn]);
            }
        }
        assert_division(level2Division, subGraphFrom0);
        auto distancesLevel2 = getDistancebetweenBoundaryNodesLevel2(
            subGraphFrom0, level2Division, boundaryNodesOfComponent);

        for (auto boundaryNode : boundaryNodesOfComponent) {
            auto shortestDistances = mainThrust(subGraphFrom0, level2Division, distancesLevel2,
                boundaryNode, boundaryNodesOfComponent);

            for (int i = 0; i < shortestDistances.size(); i++) {
                if (shortestDistances[i] != -1) {
                    result[std::make_pair(nodes[i], nodes[boundaryNode])] = shortestDistances[i];
                    result[std::make_pair(nodes[boundaryNode], nodes[i])] = shortestDistances[i];
                }
            }
        }
    }
    return result;
}

pairDistance_t getDistancebetweenBoundaryNodesLevel1(
    NetworKit::Graph& graph, nodeSubsets_t& division, int r2, int c) {
    std::vector<pairDistance_t> resultPerRegion(division.size());
    pairDistance_t result;

    std::vector<int> regionCount(graph.numberOfNodes(), 0);
    std::vector<int> isLevel1Boundary(graph.numberOfNodes(), 0);

    for (int i = 0; i < division.size(); i++) {
        for (auto node : division[i]) {
            regionCount[node]++;
        }
    }
    for (int i = 0; i < graph.numberOfNodes(); i++) {
        if (regionCount[i] > 1) {
            isLevel1Boundary[i] = true;
        }
    }

    for (int r = 0; r < division.size(); r++) {
        auto& region = division[r];
        auto subGraph = NetworKit::GraphTools::subgraphFromNodes(
            graph, std::unordered_set<NetworKit::node>{region.begin(), region.end()});
        std::vector<NetworKit::node> boundary;

        for (auto node : region) {
            if (isLevel1Boundary[node]) {
                boundary.push_back(node);
            }
        }

        resultPerRegion[r] = getDistanceInRegionLevel1(subGraph, boundary, r2, c);
    }
    for (auto pairs : resultPerRegion) {
        for (auto [pair, distance] : pairs) {
            if (result.contains(pair)) {
                result[pair] = std::min(result[pair], distance);
            } else {
                result[pair] = distance;
            }
        }
    }
    return result;
}

void FredericksonPlanarSSSP::run() {
    normal_graph = convertToMaxDegree3(graph);
    // parameters described by the paper
    // int r1 = log(normal_graph.numberOfNodes());
    // int r2 = (log(log(normal_graph.numberOfNodes())))^2

    auto divisionLevel1 = findSuitableRDivision(normal_graph, r1, c);
    assert_division(divisionLevel1, normal_graph);

    std::vector<std::vector<int>> nodeRegions(normal_graph.numberOfNodes());
    for (int i = 0; i < divisionLevel1.size(); i++) {
        for (auto node : divisionLevel1[i]) {
            nodeRegions[node].push_back(i);
        }
    }
    auto distancesLevel1 =
        getDistancebetweenBoundaryNodesLevel1(normal_graph, divisionLevel1, r2, c);

    auto shortestDistance = mainThrust(normal_graph, divisionLevel1, distancesLevel1, source, {});

    if (nodeRegions[target].size() > 1) {
        distanceToTarget = shortestDistance[target];
    } else {  // mop up is required to calculated distance to target node that is not a boundary
              // node
        distanceToTarget =
            mopUp(shortestDistance, divisionLevel1, normal_graph, nodeRegions[target][0], target);
    }

    hasRun = true;
    return;
}

} /* namespace Koala */
