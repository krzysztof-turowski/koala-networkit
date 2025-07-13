#include "planar_sssp/HenzingerPlanarSSSP.hpp"

#include <algorithm>
#include <cmath>
#include <networkit/distance/Dijkstra.hpp>
#include <networkit/distance/MultiTargetDijkstra.hpp>
#include <networkit/graph/BFS.hpp>
#include <networkit/graph/DFS.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <set>
#include <unordered_map>
#include <utility>

#include "planar_sssp/PlanarUtils.hpp"
#include "planar_sssp/SuitableRDivision.hpp"

namespace Koala {
static const int INF = std::numeric_limits<int>::max();

void HenzingerPlanarSSSP::mainThrust() {
    while (!mainQ.empty()) {
        int currentRegion = mainQ.minItem();
        std::set<int> updatedRegions;
        updatedRegions.insert(currentRegion);

        for (int times = r; times > 0; times--) {
            if (Q[currentRegion].empty()) continue;
            int distance = Q[currentRegion].minKey();
            NetworKit::node currentNode = Q[currentRegion].minItem();
            Q[currentRegion].deactivateItem(currentNode);

            for (auto [v, weight] : normal_graph.weightNeighborRange(currentNode)) {
                if (d[v] <= distance + weight) continue;
                // edge goes outside the current region and will be updated later
                auto nodeRegions = std::find(regions[v].begin(), regions[v].end(), currentRegion);
                if (nodeRegions == regions[v].end()) continue;

                d[v] = distance + weight;
                for (auto region : regions[v]) {
                    Q[region].updateKey(v, d[v]);
                    updatedRegions.insert(region);
                }
            }
        }

        updatedRegions.erase(currentRegion);
        if (!Q[currentRegion].empty()) {
            mainQ.updateKey(currentRegion, Q[currentRegion].minKey());
        } else {
            mainQ.deactivateItem(currentRegion);
        }

        for (auto region : updatedRegions) {
            mainQ.updateKey(region, Q[region].minKey());
        }
    }
}

void HenzingerPlanarSSSP::initializeQueues(nodeSubsets_t& division) {
    Q.assign(division.size(), {});
    for (int i = 0; i < division.size(); i++) {
        auto node = std::find(division[i].begin(), division[i].end(), source);
        if (node == division[i].end()) continue;
        Q[i].insert(source, 0);
        mainQ.insert(i, 0);
    }
}

void HenzingerPlanarSSSP::run() {
    normal_graph = convertToMaxDegree3(graph, true);
    auto graphForDivision = convertDirectedGraphToUndirected(normal_graph);
    r = log(graphForDivision.numberOfNodes());
    int c = 3;  // arbitrary parameter. Bounds number of boundary nodes in region of division
    // int r4 = r * r * r * r; //log(n)^4
    int r4 = 25;
    auto division = findSuitableRDivision(graphForDivision, r4, 3);
    assert_division(division, graphForDivision);
    numOfRegions = division.size();

    d.assign(normal_graph.numberOfNodes(), INF);
    d[source] = 0;

    isBoundary.assign(normal_graph.numberOfNodes(), 0);
    regions.assign(normal_graph.numberOfNodes(), {});
    for (int i = 0; i < division.size(); i++) {
        for (auto node : division[i]) {
            regions[node].push_back(i);
        }
    }
    for (int i = 0; i < normal_graph.numberOfNodes(); i++) {
        if (regions.size() > 1) {
            isBoundary[i] = 1;
        }
    }

    initializeQueues(division);

    mainThrust();

    distanceToTarget = d[target];
    hasRun = true;
    return;
}
}  // namespace Koala
