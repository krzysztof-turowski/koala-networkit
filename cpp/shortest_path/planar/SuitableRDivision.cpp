#include "shortest_path/planar/SuitableRDivision.hpp"

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

namespace Koala {
NetworKit::node root_node = 0;
std::vector<NetworKit::node> parent;
std::vector<NetworKit::count> depth;
std::vector<bool> visited;
int c;

struct Separator {
    std::unordered_set<NetworKit::node> A;
    std::unordered_set<NetworKit::node> B;
    std::unordered_set<NetworKit::node> C;
};

void inside_DFS(
    NetworKit::node u, planar_embedding_t& graph, std::unordered_set<NetworKit::node>& inside) {
    if (visited[u]) return;
    inside.insert(u);
    visited[u] = true;
    for (auto v : graph[u]) {
        inside_DFS(v, graph, inside);
    }
    return;
}

std::pair<std::unordered_set<NetworKit::node>, std::unordered_set<NetworKit::node>> count_nodes(
    NetworKit::Edge edge, planar_embedding_t& graph) {
    // TODO(289Adam289): improve cycle selection to be linear
    std::unordered_set<NetworKit::node> inside, cycle;

    for (auto& p : graph) {
        visited[p.first] = 0;
    }

    NetworKit::node last;
    NetworKit::node current = edge.u;
    NetworKit::node root_right_child;
    NetworKit::node local_root_node = NetworKit::none;

    while (current != NetworKit::none) {
        cycle.insert(current);
        visited[current] = true;
        current = parent[current];
    }
    current = edge.v;
    while (!cycle.contains(parent[current])) {
        cycle.insert(current);
        visited[current] = true;
        current = parent[current];
    }
    root_right_child = current;
    cycle.insert(current);
    visited[current] = true;
    local_root_node = parent[current];
    current = parent[local_root_node];
    while (current != NetworKit::none) {
        cycle.erase(current);
        visited[current] = false;
        current = parent[current];
    }

    // LEFT
    current = edge.u;
    last = edge.v;
    while (current != local_root_node) {
        auto start =
            std::find(graph[current].begin(), graph[current].end(), last) - graph[current].begin();
        auto end = std::find(graph[current].begin(), graph[current].end(), parent[current]) -
                   graph[current].begin();
        for (NetworKit::index i = (start + 1) % graph[current].size(); i != end;
            i = (i + 1) % graph[current].size()) {
            inside_DFS(graph[current][i], graph, inside);
        }
        last = current;
        current = parent[current];
    }

    // ROOT
    auto start = std::find(graph[local_root_node].begin(), graph[local_root_node].end(), last) -
                 graph[local_root_node].begin();
    auto end =
        std::find(graph[local_root_node].begin(), graph[local_root_node].end(), root_right_child) -
        graph[local_root_node].begin();
    for (NetworKit::index i = (start + 1) % graph[local_root_node].size(); i != end;
        i = (i + 1) % graph[local_root_node].size()) {
        inside_DFS(graph[local_root_node][i], graph, inside);
    }

    return {inside, cycle};
}

Separator find_separator(NetworKit::Graph& G) {
    root_node = *(G.nodeRange().begin());
    NetworKit::Graph maximal_graph = PlanarGraphTools::makeMaximalPlanar(G);
    NetworKit::count number_of_nodes = maximal_graph.numberOfNodes();

    for (auto node : maximal_graph.nodeRange()) {
        parent[node] = NetworKit::none;
        depth[node] = NetworKit::none;
        visited[node] = 0;
    }
    std::vector<NetworKit::Edge> tree_edges;
    std::vector<NetworKit::Edge> non_tree_edges;

    std::queue<std::pair<NetworKit::node, NetworKit::count>> queue;

    auto planar_embeding = PlanarGraphTools::findPlanarEmbeding(maximal_graph);
    parent[root_node] = NetworKit::none;
    depth[root_node] = 0;
    visited[root_node] = true;
    queue.push({root_node, 0});

    while (!queue.empty()) {
        auto [node, d] = queue.front();
        queue.pop();

        for (auto n : planar_embeding[node]) {
            if (visited[n]) {
                if (depth[node] <= depth[n]) {
                    non_tree_edges.push_back({node, n});
                }
                continue;
            }
            parent[n] = node;
            visited[n] = true;
            depth[n] = d + 1;
            tree_edges.push_back({node, n});
            queue.push({n, d + 1});
        }
    }

    std::unordered_set<NetworKit::node> inside, cycle, outside;
    for (auto [u, v] : non_tree_edges) {
        std::tie(inside, cycle) = count_nodes({u, v}, planar_embeding);
        NetworKit::count outside_count = number_of_nodes - inside.size() - cycle.size();

        assert(inside.size() + cycle.size() <= number_of_nodes);

        if (outside_count * 3 <= 2 * number_of_nodes && inside.size() * 3 <= 2 * number_of_nodes) {
            maximal_graph.forNodes([&](NetworKit::node n) {
                if (inside.find(n) == inside.end() && cycle.find(n) == cycle.end()) {
                    outside.insert(n);
                }
            });

            for (auto [a, b] : G.edgeRange()) {
                if (cycle.contains(a) || cycle.contains(b)) continue;

                // no edge that connect inside with outside
                assert((inside.contains(a) && inside.contains(b)) ||
                       (outside.contains(a) && outside.contains(b)));
            }
            assert(inside.size() + outside.size() + cycle.size() == G.numberOfNodes());

            return {inside, outside, cycle};
        }
    }

    throw std::runtime_error("Acording to Lipton and Tarjan Lemma 2 Should not have happened!!!");
}

std::vector<std::vector<NetworKit::node>> post_processing(
    Separator& separator, NetworKit::Graph& G, std::vector<int>& is_boundary) {
    auto [A, B, C] = separator;
    std::unordered_set<NetworKit::node> C_prime;
    std::unordered_set<NetworKit::node> C_bis;
    NetworKit::Graph copy_G(G);

    for (auto c : C) {
        bool is_prime = true;
        for (auto x : G.neighborRange(c)) {
            if (A.find(x) != A.end() || B.find(x) != B.end()) {
                is_prime = false;
                C_bis.insert(c);
                copy_G.removeNode(c);
                break;
            }
        }

        if (is_prime) {
            C_prime.insert(c);
        }
    }
    NetworKit::ConnectedComponents cc(copy_G);
    cc.run();
    node_subsets_t components = cc.getComponents();

    // the procedure was vaguely described in the paper
    // what's important is that the process of moving vertices from C_bis to coneccted componets is
    // greedy it's easy to prove that bellow prcess will created actuall connected subsets it
    // requires observation about structure of subgraph C_bis. C_bis is a uninion of unconnected
    // paths.

    std::unordered_map<NetworKit::node, NetworKit::count> moved_node_map;
    std::unordered_set<NetworKit::node> boundary_nodes;
    for (auto c : C_bis) {
        std::unordered_set<NetworKit::count> neighbor_components;
        for (auto x : G.neighborRange(c)) {
            if (C_bis.contains(x)) {
                if (moved_node_map.find(x) != moved_node_map.end()) {
                    neighbor_components.insert(moved_node_map[x]);
                }
                continue;
            }
            neighbor_components.insert(cc.componentOfNode(x));
        }

        if (neighbor_components.size() == 1) {
            moved_node_map[c] = *(neighbor_components.begin());
            components[moved_node_map[c]].push_back(c);
        } else {
            boundary_nodes.insert(c);
        }
    }

    std::unordered_map<NetworKit::node, int> new_boundary_node;
    for (auto node : boundary_nodes) {
        is_boundary[node] = 1;
        new_boundary_node[node] = 1;
        for (auto x : G.neighborRange(node)) {
            NetworKit::count component;

            if (C_bis.contains(x)) {
                if (moved_node_map.find(x) == moved_node_map.end()) {
                    continue;
                }
                component = moved_node_map[x];
            } else {
                component = cc.componentOfNode(x);
            }

            if (components[component].back() != node) {
                components[component].push_back(node);
            }
            moved_node_map[node] = component;
        }
    }

    // Assert will be removed later

    std::unordered_map<NetworKit::node, std::vector<NetworKit::index>> regions_of_node;
    for (NetworKit::index i = 0; i < components.size(); i++) {
        for (auto node : components[i]) {
            regions_of_node[node].push_back(i);
        }
    }

    for (auto& [k, v] : regions_of_node) {
        NetworKit::count count_neighbors = 0;
        for (auto nei : G.neighborRange(k)) {
            count_neighbors++;
        }
        assert(v.size() > 0);
        assert(v.size() <= count_neighbors);
        if (new_boundary_node[k]) assert(v.size() > 1);
    }

    for (auto [u, v] : G.edgeRange()) {
        bool has_common = std::any_of(
            regions_of_node[u].begin(), regions_of_node[u].end(), [&](NetworKit::index x) {
                return std::find(regions_of_node[v].begin(), regions_of_node[v].end(), x) !=
                       regions_of_node[v].end();
            });
        assert(has_common && "edge must be in at least one region fully");
        // in other words baundry nodes are in all regions they are connected to
    }

    // Assert ends here

    return components;
}

node_subsets_t create_connected_sets(NetworKit::Graph& subgraph, std::vector<int>& is_boundary) {
    node_subsets_t result;
    NetworKit::ConnectedComponents cc(subgraph);
    cc.run();
    node_subsets_t components = cc.getComponents();
    if (components.size() == 1) {  // subgraph is a connected componet
        auto separator = find_separator(subgraph);
        return post_processing(separator, subgraph, is_boundary);
    }

    NetworKit::count largest_component_number = NetworKit::none;
    for (NetworKit::index i = 0; i < components.size(); i++) {
        if (components[i].size() * 3 > subgraph.numberOfNodes() * 2) {
            largest_component_number = i;
        } else {
            result.push_back(std::move(components[i]));
        }
    }
    if (largest_component_number ==
        NetworKit::none) {  // all components ale smaller thatn 2/3 of number of
                            // nodes no need for using separater theorem
        return result;
    }

    // finding separator of the largest component (it size exceeds 2/3 off all nodes)
    std::unordered_set<NetworKit::node> largest_component(
        components[largest_component_number].begin(), components[largest_component_number].end());
    NetworKit::Graph connected_graph =
        NetworKit::GraphTools::subgraphFromNodes(subgraph, largest_component);
    auto separator = find_separator(connected_graph);

    auto conected_subsets = post_processing(separator, connected_graph, is_boundary);
    for (NetworKit::index i = 0; i < conected_subsets.size(); i++) {
        result.push_back(std::move(conected_subsets[i]));
    }

    return result;
}

void fix_division(node_subsets_t& division, NetworKit::Graph& graph) {
    std::vector<std::vector<NetworKit::index>> components_of_node(graph.numberOfNodes());

    for (NetworKit::index i = 0; i < division.size(); i++) {
        for (auto node : division[i]) {
            components_of_node[node].push_back(i);
        }
    }

    for (auto node : graph.nodeRange()) {
        if (components_of_node[node].size() > 3) {
            NetworKit::index component_to_remove = NetworKit::none;

            for (auto c : components_of_node[node]) {
                NetworKit::count count_neighbors = 0;
                NetworKit::node neighbor = NetworKit::none;
                for (auto nei : graph.neighborRange(node)) {
                    if (std::find(components_of_node[nei].begin(), components_of_node[nei].end(),
                            c) != components_of_node[nei].end()) {
                        count_neighbors++;
                        neighbor = nei;
                    }
                }

                if (count_neighbors == 0) {
                    component_to_remove = c;
                    break;
                }
                if (count_neighbors == 1) {
                    std::vector<NetworKit::index> rest_components;
                    for (auto cc : components_of_node[node]) {
                        if (c != cc) {
                            rest_components.push_back(cc);
                        }
                    }
                    bool has_common = std::any_of(
                        rest_components.begin(), rest_components.end(), [&](NetworKit::index x) {
                            return std::find(components_of_node[neighbor].begin(),
                                       components_of_node[neighbor].end(),
                                       x) != components_of_node[neighbor].end();
                        });
                    if (has_common) {
                        component_to_remove = c;
                        break;
                    }
                }
            }

            division[component_to_remove].erase(std::find(
                division[component_to_remove].begin(), division[component_to_remove].end(), node));
        }
    }
}

node_subsets_t make_suitable_graph_division_from_queue(
    std::queue<std::vector<NetworKit::node>>& queue, NetworKit::Graph& graph, int r,
    std::vector<int>& is_boundary) {
    int sqr = sqrt(r);
    node_subsets_t small_sets;
    node_subsets_t result;

    // dividing graph into small overlaping sets of nodes
    while (!queue.empty()) {
        auto node_set = std::move(queue.front());
        queue.pop();

        NetworKit::count number_of_boundary_nodes = 0;
        for (auto node : node_set) {
            if (is_boundary[node]) {
                number_of_boundary_nodes++;
            }
        }
        if (node_set.size() < r && number_of_boundary_nodes < c * sqr) {
            small_sets.push_back(std::move(node_set));
            continue;
        }

        std::unordered_set<NetworKit::node> nodes_for_subgraph{node_set.begin(), node_set.end()};
        auto subgraph = NetworKit::GraphTools::subgraphFromNodes(graph, nodes_for_subgraph);
        auto sets = create_connected_sets(subgraph, is_boundary);

        for (auto& s : sets) {
            number_of_boundary_nodes = 0;
            for (auto node : s) {
                if (is_boundary[node]) {
                    number_of_boundary_nodes++;
                }
            }
            if (s.size() >= r || number_of_boundary_nodes >= c * sqr) {
                queue.emplace(std::move(s));
            } else {
                small_sets.push_back(std::move(s));
            }
        }
    }

    fix_division(small_sets, graph);

    PlanarGraphTools::assertDivision(small_sets, graph);

    // merging to small components
    std::vector<std::vector<NetworKit::count>> components_per_node(graph.numberOfNodes());
    std::vector<int> is_set_used(small_sets.size(), 0);
    std::vector<NetworKit::count> number_of_boundary_nodes(small_sets.size(), 0);
    std::vector<int> used(graph.numberOfNodes(), 0);

    for (NetworKit::index i = 0; i < small_sets.size(); i++) {
        for (auto& x : small_sets[i]) {
            components_per_node[x].push_back(i);
        }
    }
    for (NetworKit::index i = 0; i < components_per_node.size(); i++) {
        if (components_per_node[i].size() > 1) {
            for (auto c : components_per_node[i]) {
                number_of_boundary_nodes[c] += 1;
            }
        }
    }

    auto is_mergable = [&](NetworKit::index c) {
        if (is_set_used[c]) return false;
        if (2 * small_sets[c].size() > r || 2 * number_of_boundary_nodes[c] > c * sqr) return false;
        return true;
    };

    auto try_merge = [&](NetworKit::index c1, NetworKit::index c2) {
        if (!is_mergable(c1)) return NetworKit::none;
        if (!is_mergable(c2)) return NetworKit::none;

        // importat to merge smaller component to the bigger
        if (small_sets[c1].size() < small_sets[c2].size()) std::swap(c1, c2);
        is_set_used[c2] = true;
        NetworKit::count good_boundary_nodes = 0;
        for (auto node : small_sets[c2]) {
            auto& node_components = components_per_node[node];
            if (std::find(node_components.begin(), node_components.end(), c1) ==
                node_components.end()) {
                small_sets[c1].push_back(node);
                if (node_components.size() > 1) {
                    good_boundary_nodes++;
                }
                std::replace(node_components.begin(), node_components.end(), c2, c1);
            } else {  // boundary between c1 and c2
                good_boundary_nodes += node_components.size() == 2 ? NetworKit::none : 1;
                node_components.erase(find(node_components.begin(), node_components.end(), c2));
            }
        }
        small_sets[c2].clear();
        number_of_boundary_nodes[c1] += good_boundary_nodes;
        assert(number_of_boundary_nodes[c1] >= 0);
        return c1;
    };

    // mergin components with common boundary nodes
    for (auto& components : components_per_node) {
        if (components.size() == 1) continue;  // interior node

        NetworKit::index c0, c1, c2;
        if (components.size() == 3) {  // boundy node of threee regions
            NetworKit::index c0 = components[0], c1 = components[1], c2 = components[2];
            try_merge(c0, c1);
            try_merge(c0, c2);
            try_merge(c1, c2);
            continue;
        }
        // boundary node of two regions
        try_merge(components[0], components[1]);
    }

    std::vector<std::vector<NetworKit::node>> regions_to_assert;
    for (auto set : small_sets) {
        if (set.size() > 0) {
            regions_to_assert.push_back(set);
        }
    }
    PlanarGraphTools::assertDivision(regions_to_assert, graph);

    NetworKit::count count_boundary_nodes = 0;
    for (auto node : graph.nodeRange()) {
        assert(components_per_node[node].size() > 0);
        if (components_per_node[node].size() > 1) {
            count_boundary_nodes++;
        }
    }

    std::vector<std::tuple<NetworKit::index, NetworKit::index, NetworKit::index>> regions_to_merge;

    for (NetworKit::index i = 0; i < small_sets.size(); i++) {
        if (small_sets[i].size() == 0) continue;
        if (!is_mergable(i)) continue;
        std::unordered_set<NetworKit::index> region_neighbors;
        for (auto node : small_sets[i]) {
            for (auto c : components_per_node[node]) {
                if (c != i) {
                    region_neighbors.insert(c);
                }
            }
        }
        if (region_neighbors.size() > 2) {
            // regions with more than two neighbors can be moved straight to the result
            // we can delay that to simplify the logic as the sets are still good to use
            continue;
        }

        NetworKit::index first, second;
        auto it = region_neighbors.begin();
        first = *it;
        it++;
        second = it == region_neighbors.end() ? NetworKit::none : *it;

        regions_to_merge.push_back({first, second, i});
    }

    sort(regions_to_merge.begin(), regions_to_merge.end());
    // TODO(289Adam289): swap std::sort for an linear or nsqrt(logn) sort

    // greedy merging of regions with the same sets of either one or two regions.

    NetworKit::index last_component = NetworKit::none;
    std::pair<NetworKit::index, NetworKit::index> last_pair = {NetworKit::none, NetworKit::none};
    for (auto [f, s, i] : regions_to_merge) {
        if (!is_mergable(i)) continue;
        if (last_component = NetworKit::none) {
            last_component = i;
            last_pair = {f, s};
            continue;
        }
        if (std::make_pair(f, s) != last_pair) {
            last_component = i;
            last_pair = {f, s};
            continue;
        }

        NetworKit::index merged = try_merge(last_component, i);
        if (!is_mergable(merged)) {
            last_component = NetworKit::none;
            continue;
        }
        if (merged == last_component) {
            continue;
        }
        last_component = merged;
        last_pair = {f, s};
    }

    for (NetworKit::index i = 0; i < small_sets.size(); i++) {
        if (is_set_used[i]) continue;
        result.push_back(std::move(small_sets[i]));
    }

    PlanarGraphTools::assertDivision(result, graph);
    return result;
}

// we need two different starting points for algorithm
node_subsets_t make_suitable_graph_division(NetworKit::Graph& graph, int r) {
    std::vector<int> is_boundary(graph.numberOfNodes(), 0);
    std::queue<std::vector<NetworKit::node>> queue;
    auto node_range = graph.nodeRange();
    queue.push(std::vector<NetworKit::node>{node_range.begin(), node_range.end()});

    return make_suitable_graph_division_from_queue(queue, graph, r, is_boundary);
}

// FIND_CLUSTERS from [F1]
node_subsets_t find_cluster_result;

std::vector<NetworKit::node> csearch(
    NetworKit::node v, NetworKit::count z, NetworKit::Graph& graph) {
    if (visited[v]) {
        return {};
    }
    visited[v] = 1;
    std::vector<NetworKit::node> clust(1, v);

    for (auto neigh : graph.neighborRange(v)) {
        auto neighbor_cluster = csearch(neigh, z, graph);
        for (auto n : neighbor_cluster) {
            clust.push_back(n);
        }
    }
    if (clust.size() < z) {
        return clust;
    } else {
        find_cluster_result.push_back(clust);
        return {};
    }
}

node_subsets_t find_clusters(NetworKit::Graph& graph, NetworKit::count z) {
    find_cluster_result.clear();
    for (auto node : graph.nodeRange()) {
        visited[node] = 0;
    }

    auto csearch_result = csearch(*graph.nodeRange().begin(), z, graph);

    for (auto node : csearch_result) {
        find_cluster_result.back().push_back(node);
    }

    return find_cluster_result;
}

node_subsets_t get_division_from_clusters(
    node_subsets_t& clusters, NetworKit::Graph& graph, int r) {
    std::vector<NetworKit::count> cluster_numbers(graph.numberOfNodes());
    for (NetworKit::index i = 0; i < clusters.size(); i++) {
        for (auto node : clusters[i]) {
            cluster_numbers[node] = i;
        }
    }
    // build new shrinked graph based on clusters
    NetworKit::Graph shrinked_graph(clusters.size());

    graph.forNodes([&](auto v) {
        for (auto u : graph.neighborRange(v)) {
            if (cluster_numbers[v] != cluster_numbers[u]) {
                shrinked_graph.addEdge(cluster_numbers[v], cluster_numbers[u]);
            }
        }
    });
    shrinked_graph.removeMultiEdges();

    int sqr = sqrt(r);
    node_subsets_t division;

    std::vector<int> is_boundary(shrinked_graph.numberOfNodes(), 0);
    std::queue<std::vector<NetworKit::node>> queue;
    auto node_range = shrinked_graph.nodeRange();
    queue.push(std::vector<NetworKit::node>{node_range.begin(), node_range.end()});

    // dividing graph into small overlaping sets of nodes
    while (!queue.empty()) {
        auto node_set = std::move(queue.front());
        queue.pop();

        NetworKit::count number_of_boundary_nodes = 0;
        for (auto node : node_set) {
            if (is_boundary[node]) {
                number_of_boundary_nodes++;
            }
        }
        if (node_set.size() < r && number_of_boundary_nodes < c * sqr) {
            division.push_back(std::move(node_set));
            continue;
        }

        std::unordered_set<NetworKit::node> nodes_for_subgraph{node_set.begin(), node_set.end()};
        auto subgraph =
            NetworKit::GraphTools::subgraphFromNodes(shrinked_graph, nodes_for_subgraph);
        auto sets = create_connected_sets(subgraph, is_boundary);

        for (auto& s : sets) {
            number_of_boundary_nodes = 0;
            for (auto node : s) {
                if (is_boundary[node]) {
                    number_of_boundary_nodes++;
                }
            }
            if (s.size() >= r || number_of_boundary_nodes >= c * sqr) {
                queue.emplace(std::move(s));
            } else {
                division.push_back(std::move(s));
            }
        }
    }

    std::vector<NetworKit::index> region_of_node(shrinked_graph.numberOfNodes(), NetworKit::none);
    std::vector<NetworKit::count> count_regions_of_node(shrinked_graph.numberOfNodes(), 0);

    for (NetworKit::index i = 0; i < division.size(); i++) {
        for (auto& node : division[i]) {
            count_regions_of_node[node] += 1;
            region_of_node[node] = i;
        }
    }

    node_subsets_t result_with_zeros(division.size(), std::vector<NetworKit::node>{});

    // Unshrink the graph
    for (NetworKit::index i = 0; i < count_regions_of_node.size(); i++) {
        assert(region_of_node[i] != NetworKit::none);
        // boundary vertex
        if (count_regions_of_node[i] >= 2) {
            result_with_zeros.push_back({});
            for (auto node : clusters[i]) {
                result_with_zeros.back().push_back(node);
            }
        } else {  // interior
            for (auto node : clusters[i]) {
                result_with_zeros[region_of_node[i]].push_back(node);
            }
        }
    }

    node_subsets_t result;

    for (NetworKit::index i = 0; i < result_with_zeros.size(); i++) {
        if (result_with_zeros[i].size() > 0) {
            result.push_back(std::move(result_with_zeros[i]));
        }
    }

    return result;
}

node_subsets_t fix_graph_division(
    node_subsets_t& division, NetworKit::Graph& graph, std::vector<int>& is_boundary) {
    // The procedure is not described in the paper "Infer boundary vertices and slightly
    // expanded(sic!) regions that share these vertices". So I guess we can do it in a greedy
    // fashion.

    std::vector<std::array<NetworKit::index, 3>> region(
        graph.numberOfNodes(), {NetworKit::none, NetworKit::none, NetworKit::none});
    for (NetworKit::index i = 0; i < division.size(); i++) {
        for (auto& node : division[i]) {
            region[node][0] = i;
        }
    }

    graph.forEdges([&](NetworKit::node u, NetworKit::node v) {
        if (is_boundary[u] && is_boundary[v]) {
            return;
        }
        if (is_boundary[u] || is_boundary[v]) {
            if (is_boundary[v]) std::swap(u, v);
            for (auto r : region[u]) {
                if (r == region[v][0]) return;
            }
            region[u][2] = region[v][0];
            return;
        }
        if (region[u][0] == region[v][0]) return;

        is_boundary[u] = 1;
        region[u][1] = region[v][0];
    });

    graph.forEdges([&](NetworKit::node u, NetworKit::node v) {
        if (is_boundary[u] && is_boundary[v]) {
            bool has_common = false;

            for (auto cu : region[u]) {
                if (cu == NetworKit::none) break;
                for (auto cv : region[v]) {
                    if (cv == cu) has_common = true;
                }
            }
            if (!has_common) {
                if (region[u][2] == NetworKit::none) {
                    region[u][2] = region[v][0];
                } else {
                    region[v][2] = region[u][0];
                }
            }
        }
    });

    node_subsets_t result(division.size());

    for (NetworKit::index i = 0; i < region.size(); i++) {
        for (auto reg : region[i]) {
            if (reg == NetworKit::none) continue;
            result[reg].push_back(i);
        }
    }

    return result;
}

// create suitable r-division quickly.
node_subsets_t findSuitableRDivision(NetworKit::Graph& graph, int r, int const_C) {
    // values used globally in dfs
    parent.resize(graph.numberOfNodes());
    depth.resize(graph.numberOfNodes());
    visited.resize(graph.numberOfNodes());
    c = const_C;
    int sqrt_r = sqrt(r);

    auto clusters = find_clusters(graph, sqrt_r);

    node_subsets_t division = get_division_from_clusters(clusters, graph, r);

    for (auto& region : division) {
        assert(region.size() > 0);
    }

    std::vector<int> is_boundary(graph.numberOfNodes(), 0);
    division = fix_graph_division(division, graph, is_boundary);

    PlanarGraphTools::assertDivision(division, graph);

    std::queue<std::vector<NetworKit::node>> queue;
    for (auto& region : division) {
        queue.push(std::move(region));
    }
    return make_suitable_graph_division_from_queue(queue, graph, r, is_boundary);
}
}  // namespace Koala
