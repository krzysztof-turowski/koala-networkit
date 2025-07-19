#include "shortest_path/HenzingerPlanarSSSP.hpp"

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

#include "graph/GraphTools.hpp"
#include "graph/PlanarGraphTools.hpp"
#include "shortest_path/planar/SuitableRDivision.hpp"

namespace Koala {
static const NetworKit::edgeweight INF = std::numeric_limits<NetworKit::edgeweight>::max();

void HenzingerPlanarSSSP::main_thrust() {
    while (!main_Q.empty()) {
        NetworKit::index current_region = main_Q.minimum_item();
        std::set<NetworKit::index> updated_regions;
        updated_regions.insert(current_region);

        for (NetworKit::index times = r; times > 0; times--) {
            if (Q[current_region].empty()) continue;
            NetworKit::edgeweight distance = Q[current_region].minimum_key();
            NetworKit::node current_node = Q[current_region].minimum_item();
            Q[current_region].deactivate(current_node);

            for (auto [v, weight] : normal_graph.weightNeighborRange(current_node)) {
                if (d[v] <= distance + weight) continue;
                // edge goes outside the current region and will be updated later
                auto node_regions = std::find(regions[v].begin(), regions[v].end(), current_region);
                if (node_regions == regions[v].end()) continue;

                d[v] = distance + weight;
                for (auto region : regions[v]) {
                    Q[region].update(v, d[v]);
                    updated_regions.insert(region);
                }
            }
        }

        updated_regions.erase(current_region);
        if (!Q[current_region].empty()) {
            main_Q.update(current_region, Q[current_region].minimum_key());
        } else {
            main_Q.deactivate(current_region);
        }

        for (auto region : updated_regions) {
            main_Q.update(region, Q[region].minimum_key());
        }
    }
}

void HenzingerPlanarSSSP::initialize_queues(node_subsets_t& division) {
    Q.assign(division.size(), {});
    for (NetworKit::index i = 0; i < division.size(); i++) {
        auto node = std::find(division[i].begin(), division[i].end(), source);
        if (node == division[i].end()) continue;
        Q[i].push(source, 0);
        main_Q.push(i, 0);
    }
}

void HenzingerPlanarSSSP::run() {
    normal_graph = PlanarGraphTools::convertToMaxDegree3(graph, true);
    auto graph_for_division = GraphTools::convertDirectedGraphToUndirected(normal_graph);
    r = log(graph_for_division.numberOfNodes());
    int c = 3;  // arbitrary parameter. Bounds number of boundary nodes in region of division
    // int r4 = r * r * r * r; //log(n)^4
    int r4 = 25;
    auto division = findSuitableRDivision(graph_for_division, r4, 3);
    PlanarGraphTools::assertDivision(division, graph_for_division);
    number_of_regions = division.size();

    d.assign(normal_graph.numberOfNodes(), INF);
    d[source] = 0;

    is_boundary.assign(normal_graph.numberOfNodes(), 0);
    regions.assign(normal_graph.numberOfNodes(), {});
    for (NetworKit::index i = 0; i < division.size(); i++) {
        for (auto node : division[i]) {
            regions[node].push_back(i);
        }
    }
    for (NetworKit::index i = 0; i < normal_graph.numberOfNodes(); i++) {
        if (regions.size() > 1) {
            is_boundary[i] = 1;
        }
    }
    initialize_queues(division);

    main_thrust();

    for (auto node : graph.nodeRange()) {
        distances[node] = d[node];
    }

    hasRun = true;
    return;
}
}  // namespace Koala
