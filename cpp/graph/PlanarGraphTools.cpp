#include "graph/PlanarGraphTools.hpp"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/biconnected_components.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>
#include <boost/graph/make_biconnected_planar.hpp>
#include <boost/graph/make_maximal_planar.hpp>
#include <networkit/graph/BFS.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace boost;

// Boost graph typedefs
using boost_graph_t = adjacency_list<vecS, vecS, undirectedS,
    property<vertex_index_t, NetworKit::edgeweight>, property<edge_index_t, NetworKit::edgeweight>>;

// Boost planar embedding typedefs based on
// https://www.boost.org/doc/libs/1_53_0/libs/graph/example/straight_line_drawing.cpp
using embedding_storage_t = std::vector<std::vector<graph_traits<boost_graph_t>::edge_descriptor>>;
using embedding_t = boost::iterator_property_map<embedding_storage_t::iterator,
    property_map<boost_graph_t, vertex_index_t>::type>;

namespace Koala {

namespace PlanarGraphTools {

std::pair<boost_graph_t, std::vector<NetworKit::node>> convert_networKit_to_boost(
    const NetworKit::Graph& networkit_graph) {
    std::vector<NetworKit::node> node_map;
    std::unordered_map<NetworKit::node, NetworKit::node> reverse_map;
    networkit_graph.forNodes([&](NetworKit::node u) {
        node_map.push_back(u);
        reverse_map[u] = node_map.size() - 1;
    });
    boost_graph_t boost_graph(networkit_graph.numberOfNodes());
    networkit_graph.forEdges([&](NetworKit::node u, NetworKit::node v) {
        add_edge(reverse_map[u], reverse_map[v], boost_graph);
    });
    return {boost_graph, node_map};
}

NetworKit::Graph convert_boost_to_networKit(
    boost_graph_t G, std::vector<NetworKit::node> node_map, NetworKit::Graph graph) {
    NetworKit::Graph result(graph);

    for (auto [ei, ei_end] = edges(G); ei != ei_end; ++ei) {
        auto u = source(*ei, G);
        auto v = target(*ei, G);

        result.addEdge(node_map[u], node_map[v]);
    }
    return result;
}

planar_embedding_t findPlanarEmbedding(const NetworKit::Graph& G, bool verbose = false) {
    auto [boost_graph, node_map] = convert_networKit_to_boost(G);
    embedding_storage_t embedding_storage(num_vertices(boost_graph));
    embedding_t embedding(embedding_storage.begin(), get(vertex_index, boost_graph));

    if (!boyer_myrvold_planarity_test(
            boyer_myrvold_params::graph = boost_graph,
            boyer_myrvold_params::embedding = embedding)) {
        throw std::runtime_error("Graph have to be planar!");
    }

    planar_embedding_t result;
    for (NetworKit::index i = 0; i < embedding_storage.size(); i++) {
        for (auto& edge : embedding_storage[i]) {
            NetworKit::node node = (
                source(edge, boost_graph) == i
                ? target(edge, boost_graph) : source(edge, boost_graph));
            if (result[node_map[i]].empty() || result[node_map[i]].back() != node_map[node]) {
                result[node_map[i]].push_back(node_map[node]);
            }
        }
    }
    return result;
}

NetworKit::Graph makeMaximalPlanar(NetworKit::Graph& G) {
    typedef adjacency_list<vecS, vecS, undirectedS, property<vertex_index_t, NetworKit::edgeweight>,
        property<edge_index_t, NetworKit::edgeweight>>
        graph;
    auto [g, node_map] = convert_networKit_to_boost(G);

    property_map<graph, edge_index_t>::type e_index = get(edge_index, g);
    graph_traits<graph>::edges_size_type edge_count = 0;
    graph_traits<graph>::edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) put(e_index, *ei, edge_count++);

    typedef std::vector<graph_traits<graph>::edge_descriptor> vec_t;
    std::vector<vec_t> embedding(num_vertices(g));
    boyer_myrvold_planarity_test(
        boyer_myrvold_params::graph = g, boyer_myrvold_params::embedding = &embedding[0]);

    make_biconnected_planar(g, &embedding[0]);

    edge_count = 0;
    for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) put(e_index, *ei, edge_count++);

    boyer_myrvold_planarity_test(
        boyer_myrvold_params::graph = g, boyer_myrvold_params::embedding = &embedding[0]);

    make_maximal_planar(g, &embedding[0]);

    // Re-initialize the edge index, since we just added a few edges
    edge_count = 0;
    for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) put(e_index, *ei, edge_count++);

    return convert_boost_to_networKit(g, node_map, G);
}

NetworKit::Graph convertToMaxDegree3(NetworKit::Graph& G, bool directed = false) {
    planar_embedding_t embedding = findPlanarEmbedding(G);
    NetworKit::count result_size = G.numberOfNodes();

    for (auto [i, node_embeding] : embedding) {
        if (node_embeding.size() > 3) {
            result_size += node_embeding.size() - 1;
        }
    }

    if (result_size == G.numberOfNodes()) {
        return G;
    }

    std::vector<std::unordered_map<size_t, size_t>> vertex_maps(
        G.numberOfNodes(), std::unordered_map<size_t, size_t>());

    size_t empty_vertex = G.numberOfNodes();
    NetworKit::Graph result(result_size, true, directed);

    for (size_t i = 0; i < embedding.size(); ++i) {
        if (embedding[i].size() > 3) {
            size_t t = embedding[i][0];
            vertex_maps[i][t] = i;
            vertex_maps[i][empty_vertex] = 0;

            result.addEdge(i, empty_vertex, 0);
            for (NetworKit::index j = 1; j < embedding[i].size(); ++j) {
                t = embedding[i][j];
                vertex_maps[i][t] = empty_vertex;
                if (j < embedding[i].size() - 1) {
                    result.addEdge(empty_vertex, empty_vertex + 1, 0);
                }
                empty_vertex++;
            }

            result.addEdge(empty_vertex - 1, i, 0);
        }
    }

    for (const auto& edge : G.edgeWeightRange()) {
        NetworKit::node u = edge.u;
        NetworKit::node v = edge.v;
        NetworKit::edgeweight ew = edge.weight;

        size_t s = v;
        size_t t = u;
        if (vertex_maps[u].find(v) != vertex_maps[u].end()) {
            t = vertex_maps[u][v];
        }
        if (vertex_maps[v].find(u) != vertex_maps[v].end()) {
            s = vertex_maps[v][u];
        }
        result.addEdge(t, s, ew);
    }

    for (auto node : result.nodeRange()) {
        NetworKit::count count = 0;
        for (auto neighbor : result.neighborRange(node)) {
            count++;
        }
        assert(count < 4);
    }

    return result;
}

void assertDivision(const node_subsets_t& division, NetworKit::Graph& graph) {
    std::vector<std::vector<NetworKit::index>> components_of_node(graph.numberOfNodes());

    for (NetworKit::index i = 0; i < division.size(); i++) {
        for (auto node : division[i]) {
            components_of_node[node].push_back(i);
        }
    }

    for (auto& components : components_of_node) {
        assert(components.size() > 0);
        assert(components.size() < 4);
    }

    for (const auto& [u, v] : graph.edgeRange()) {
        bool has_common = std::any_of(
            components_of_node[u].begin(), components_of_node[u].end(), [&](NetworKit::index x) {
                return std::find(components_of_node[v].begin(), components_of_node[v].end(), x) !=
                       components_of_node[v].end();
            });
        assert(has_common && "edge must be in at least one region fully");
    }
}

}  // namespace PlanarGraphTools

}  // namespace Koala
