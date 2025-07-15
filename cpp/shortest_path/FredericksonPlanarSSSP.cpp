#include "shortest_path/FredericksonPlanarSSSP.hpp"

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

#include "graph/PlanarGraphTools.hpp"
#include "shortest_path/planar/SuitableRDivision.hpp"
#include "shortest_path/planar/TopologyBasedHeap.hpp"

namespace Koala {
pair_distance_t get_distance_between_boundary_nodes_level_2(NetworKit::Graph& graph,
    node_subsets_t& division, const std::unordered_set<NetworKit::node>& extra_boundary_nodes) {
    std::vector<pair_distance_t> result_per_region(division.size());
    pair_distance_t result;

    std::vector<int> region_count(graph.numberOfNodes(), 0);
    std::vector<int> is_boundary(graph.numberOfNodes(), 0);

    for (int i = 0; i < division.size(); i++) {
        for (auto node : division[i]) {
            region_count[node]++;
        }
    }
    for (int i = 0; i < graph.numberOfNodes(); i++) {
        if (region_count[i] > 1) {
            is_boundary[i] = true;
        }
    }

    for (int r = 0; r < division.size(); r++) {
        auto& region = division[r];
        auto subgraph = NetworKit::GraphTools::subgraphFromNodes(
            graph, std::unordered_set<NetworKit::node>{region.begin(), region.end()});
        std::vector<NetworKit::node> boundary;

        for (auto node : region) {
            if (is_boundary[node] ||
                (extra_boundary_nodes.find(node) != extra_boundary_nodes.end())) {
                boundary.push_back(node);
            }
        }

        for (int i = 0; i < boundary.size(); i++) {
            NetworKit::MultiTargetDijkstra dijkstra(
                subgraph, boundary[i], boundary.begin(), boundary.end());
            dijkstra.run();
            auto map = dijkstra.getTargetIndexMap();
            auto distances = dijkstra.getDistances();
            for (int j = 0; j < boundary.size(); j++) {
                if (distances[map[boundary[j]]] < 0 ||
                    distances[map[boundary[j]]] > std::numeric_limits<int>::max())
                    continue;
                result_per_region[r][std::make_pair(boundary[i], boundary[j])] =
                    distances[map[boundary[j]]];
                result_per_region[r][std::make_pair(boundary[j], boundary[i])] =
                    distances[map[boundary[j]]];
            }
        }
    }

    for (auto pairs : result_per_region) {
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

void mop_up(std::vector<int>& shortest_distance, node_subsets_t& regions, NetworKit::Graph& G) {
    for (auto& region : regions) {
        std::unordered_set<NetworKit::node> subgraph_nodes(region.begin(), region.end());
        NetworKit::Graph subgraph = NetworKit::GraphTools::subgraphFromNodes(G, subgraph_nodes);
        auto root = subgraph.addNode();
        for (auto node : region) {
            if (shortest_distance[node] != -1) {
                subgraph.addEdge(root, node, shortest_distance[node]);
            }
        }
        NetworKit::MultiTargetDijkstra dijkstra(subgraph, root, region.begin(), region.end());
        dijkstra.run();
        auto map = dijkstra.getTargetIndexMap();
        auto distances = dijkstra.getDistances();

        for (auto node : region) {
            if (shortest_distance[node] != -1) continue;
            shortest_distance[node] = distances[map[node]];
        }
    }
}

std::vector<int> main_thrust(NetworKit::Graph& graph, node_subsets_t& regions,
    pair_distance_t& distances, int s,
    const std::unordered_set<NetworKit::node>& extra_boudary_nodes) {
    std::vector<int> shortest_distance(graph.numberOfNodes(), -1);

    TopologyHeap heap(graph, regions, distances, s, extra_boudary_nodes);

    while (!heap.empty()) {
        auto [distance, node] = heap.top();
        if (node == -1) break;
        shortest_distance[node] = distance;
        heap.close(node);
    }

    return shortest_distance;
}

pair_distance_t get_distance_in_region_level_1(
    NetworKit::Graph& graph, std::vector<NetworKit::node>& level_1_boundary_nodes, int r2, int c) {
    pair_distance_t result;

    NetworKit::ConnectedComponents cc(graph);
    cc.run();
    node_subsets_t components = cc.getComponents();

    for (auto& component : components) {
        auto connected_subgraph =
            NetworKit::GraphTools::subgraphFromNodes(graph, {component.begin(), component.end()});
        std::vector<NetworKit::node> nodes;
        std::unordered_map<NetworKit::node, NetworKit::node> map_of_nodes;

        for (auto node : connected_subgraph.nodeRange()) nodes.push_back(node);
        for (int i = 0; i < nodes.size(); i++) map_of_nodes[nodes[i]] = i;

        NetworKit::Graph subgraph_from_0(nodes.size(), true);
        for (auto e : connected_subgraph.edgeWeightRange()) {
            subgraph_from_0.setWeight(map_of_nodes[e.u], map_of_nodes[e.v], e.weight);
        }

        node_subsets_t level_2_division;
        if (subgraph_from_0.numberOfNodes() <= r2) {
            level_2_division.push_back(std::vector<NetworKit::node>{});
            for (auto node : subgraph_from_0.nodeRange()) {
                level_2_division[0].push_back(node);
            }
        } else {
            level_2_division = findSuitableRDivision(subgraph_from_0, r2, c);
        }

        std::unordered_set<NetworKit::node> boundary_nodes_of_component;
        for (auto bn : level_1_boundary_nodes) {
            if (map_of_nodes.find(bn) != map_of_nodes.end()) {
                boundary_nodes_of_component.insert(map_of_nodes[bn]);
            }
        }
        PlanarGraphTools::assertDivision(level_2_division, subgraph_from_0);
        auto distancesLevel2 = get_distance_between_boundary_nodes_level_2(
            subgraph_from_0, level_2_division, boundary_nodes_of_component);

        for (auto boundaryNode : boundary_nodes_of_component) {
            auto shortest_distances = main_thrust(subgraph_from_0, level_2_division,
                distancesLevel2, boundaryNode, boundary_nodes_of_component);

            for (int i = 0; i < shortest_distances.size(); i++) {
                if (shortest_distances[i] != -1) {
                    result[std::make_pair(nodes[i], nodes[boundaryNode])] = shortest_distances[i];
                    result[std::make_pair(nodes[boundaryNode], nodes[i])] = shortest_distances[i];
                }
            }
        }
    }
    return result;
}

pair_distance_t get_distance_between_boundary_nodes_Level_1(
    NetworKit::Graph& graph, node_subsets_t& division, int r2, int c) {
    std::vector<pair_distance_t> result_per_region(division.size());
    pair_distance_t result;

    std::vector<int> region_count(graph.numberOfNodes(), 0);
    std::vector<int> is_level_1_boundary(graph.numberOfNodes(), 0);

    for (int i = 0; i < division.size(); i++) {
        for (auto node : division[i]) {
            region_count[node]++;
        }
    }
    for (int i = 0; i < graph.numberOfNodes(); i++) {
        if (region_count[i] > 1) {
            is_level_1_boundary[i] = true;
        }
    }

    for (int r = 0; r < division.size(); r++) {
        auto& region = division[r];
        auto subgraph = NetworKit::GraphTools::subgraphFromNodes(
            graph, std::unordered_set<NetworKit::node>{region.begin(), region.end()});
        std::vector<NetworKit::node> boundary;

        for (auto node : region) {
            if (is_level_1_boundary[node]) {
                boundary.push_back(node);
            }
        }

        result_per_region[r] = get_distance_in_region_level_1(subgraph, boundary, r2, c);
    }
    for (auto pairs : result_per_region) {
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
    normal_graph = PlanarGraphTools::convertToMaxDegree3(graph);
    // parameters described by the paper
    // int r1 = log(normal_graph.numberOfNodes());
    // int r2 = (log(log(normal_graph.numberOfNodes())))^2

    auto division_level_1 = findSuitableRDivision(normal_graph, r1, c);
    PlanarGraphTools::assertDivision(division_level_1, normal_graph);

    std::vector<std::vector<int>> node_regions(normal_graph.numberOfNodes());
    for (int i = 0; i < division_level_1.size(); i++) {
        for (auto node : division_level_1[i]) {
            node_regions[node].push_back(i);
        }
    }
    auto distances_level_1 =
        get_distance_between_boundary_nodes_Level_1(normal_graph, division_level_1, r2, c);

    auto shortest_distances =
        main_thrust(normal_graph, division_level_1, distances_level_1, source, {});

    mop_up(shortest_distances, division_level_1, normal_graph);
    for (auto node : graph.nodeRange()) {
        distances[node] = shortest_distances[node];
    }

    hasRun = true;
    return;
}

} /* namespace Koala */
