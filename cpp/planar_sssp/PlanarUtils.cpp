#include <planar_sssp/PlanarUtils.hpp> 

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>
#include <boost/graph/biconnected_components.hpp>
#include <boost/graph/make_biconnected_planar.hpp>
#include <boost/graph/make_maximal_planar.hpp>
#include <set>
#include <vector>
#include <unordered_map>
#include <utility>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/graph/BFS.hpp>

using namespace boost;

// Boost graph typedefs
using BoostGraph =
adjacency_list<vecS, vecS, undirectedS, property<vertex_index_t, int>,
    property< edge_index_t, int >>;

// Boost planar embedding typedefs based on
// https://www.boost.org/doc/libs/1_53_0/libs/graph/example/straight_line_drawing.cpp
using embedding_storage_t =
std::vector<std::vector<graph_traits<BoostGraph>::edge_descriptor>>;
using embedding_t = boost::iterator_property_map<
    embedding_storage_t::iterator,
    property_map<BoostGraph, vertex_index_t>::type>;

using pairDistance_t = std::unordered_map<std::pair<int, int>, int, boost::hash<std::pair<int, int>>>;
using nodeSubsets_t = std::vector<std::vector<NetworKit::node>>;

namespace Koala {

    std::pair<BoostGraph, std::vector<NetworKit::node>> convertNetworKitToBoost(const NetworKit::Graph& netGraph) {
        std::vector<NetworKit::node> nodeMap;
        std::unordered_map<NetworKit::node, NetworKit::node> reverseMap;
        netGraph.forNodes([&](NetworKit::node u) {
            nodeMap.push_back(u);
            reverseMap[u] = nodeMap.size() - 1;
            });
        BoostGraph boostGraph(netGraph.numberOfNodes());
        netGraph.forEdges([&](NetworKit::node u, NetworKit::node v) { add_edge(reverseMap[u], reverseMap[v], boostGraph); });
        return { boostGraph, nodeMap };
    }


    NetworKit::Graph convertBoostToNetworKit(BoostGraph G, std::vector<NetworKit::node> nodeMap, NetworKit::Graph Graph) {
        NetworKit::Graph result(Graph);

        for (auto [ei, ei_end] = edges(G); ei != ei_end; ++ei) {
            auto u = source(*ei, G);
            auto v = target(*ei, G);

            result.addEdge(nodeMap[u], nodeMap[v]);
        }
        return result;
    }

    void printEmbeding(const planar_embedding_t& embedding_storage,
        const BoostGraph& boostGraph) {
        for (auto [k, v] : embedding_storage) {
            std::cout << "Vertex " << k << ": ";
            for (auto x : v) {
                std::cout << x << " ";
            }
            std::cout << std::endl;
        }
    }

    planar_embedding_t findPlanarEmbeding(const NetworKit::Graph& G, bool verbose = false) {
        auto [boostGraph, nodeMap] = convertNetworKitToBoost(G);
        embedding_storage_t embedding_storage(num_vertices(boostGraph));
        embedding_t embedding(embedding_storage.begin(),
            get(vertex_index, boostGraph));

        if (!boyer_myrvold_planarity_test(
            boyer_myrvold_params::graph = boostGraph,
            boyer_myrvold_params::embedding = embedding)) {
            throw std::runtime_error("Graph have to be planar!");
            return;
        }

        planar_embedding_t result;

        for (int i = 0; i < embedding_storage.size(); i++) {
            for (auto& edge : embedding_storage[i]) {
                NetworKit::node node = source(edge, boostGraph) == i ? target(edge, boostGraph) : source(edge, boostGraph);
                if (result[nodeMap[i]].empty() || result[nodeMap[i]].back() != nodeMap[node]) {
                    result[nodeMap[i]].push_back(nodeMap[node]);
                }

            }
        }

        if (verbose) {
            printEmbeding(result, boostGraph);
        }

        return result;
    }

    NetworKit::Graph makeMaximalPlanar(NetworKit::Graph& G) {
        typedef adjacency_list< vecS, vecS, undirectedS,
            property< vertex_index_t, int >, property< edge_index_t, int > >
            graph;
        auto [g, nodeMap] = convertNetworKitToBoost(G);

        property_map< graph, edge_index_t >::type e_index = get(edge_index, g);
        graph_traits< graph >::edges_size_type edge_count = 0;
        graph_traits< graph >::edge_iterator ei, ei_end;
        for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
            put(e_index, *ei, edge_count++);

        typedef std::vector< graph_traits< graph >::edge_descriptor > vec_t;
        std::vector< vec_t > embedding(num_vertices(g));
        boyer_myrvold_planarity_test(boyer_myrvold_params::graph = g,
            boyer_myrvold_params::embedding = &embedding[0]);


        make_biconnected_planar(g, &embedding[0]);

        edge_count = 0;
        for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
            put(e_index, *ei, edge_count++);

        boyer_myrvold_planarity_test(boyer_myrvold_params::graph = g,
            boyer_myrvold_params::embedding = &embedding[0]);

        make_maximal_planar(g, &embedding[0]);

        // Re-initialize the edge index, since we just added a few edges
        edge_count = 0;
        for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
            put(e_index, *ei, edge_count++);

        return convertBoostToNetworKit(g, nodeMap, G);
    }

    NetworKit::Graph convertToMaxDegree3(NetworKit::Graph& G) {
        planar_embedding_t embedding = findPlanarEmbeding(G);
        NetworKit::count result_size = G.numberOfNodes();

        for (auto [i, node_embeding] : embedding) {
            if (node_embeding.size() > 3) {
                result_size += node_embeding.size() - 1;
            }
        }

        if (result_size == G.numberOfNodes()) {
            return G;
        }

        std::vector<std::unordered_map<size_t, size_t>> vertexMaps(
            G.numberOfNodes(), std::unordered_map<size_t, size_t>());

        size_t emptyVertex = G.numberOfNodes();
        NetworKit::Graph result(result_size, true, false);

        for (size_t i = 0; i < embedding.size(); ++i) {
            if (embedding[i].size() > 3) {
                size_t t = embedding[i][0];
                vertexMaps[i][t] = i;
                vertexMaps[i][emptyVertex] = 0;

                result.addEdge(i, emptyVertex, 0);
                for (int j = 1; j < embedding[i].size(); ++j) {
                    t = embedding[i][j];
                    vertexMaps[i][t] = emptyVertex;
                    if (j < embedding[i].size() - 1) {
                        result.addEdge(emptyVertex, emptyVertex + 1, 0);
                    }
                    emptyVertex++;
                }

                result.addEdge(emptyVertex - 1, i, 0);
            }
        }

        for (const auto& edge : G.edgeWeightRange()) {
            NetworKit::node u = edge.u;
            NetworKit::node v = edge.v;
            NetworKit::edgeweight ew = edge.weight;

            size_t s = v;
            size_t t = u;
            if (vertexMaps[u].find(v) != vertexMaps[u].end()) {
                t = vertexMaps[u][v];
            }
            if (vertexMaps[v].find(u) != vertexMaps[v].end()) {
                s = vertexMaps[v][u];
            }
            result.addEdge(s, t, ew);
        }

        for (auto node : result.nodeRange()) {
            int count = 0;
            for (auto nei : result.neighborRange(node)) {
                count++;
            }
            assert(count < 4);
        }

        return result;
    }

    void assert_division(const nodeSubsets_t& division, NetworKit::Graph& Graph) {
        std::vector<std::vector<int>> componentsOfNode(Graph.numberOfNodes());

        for (int i = 0; i < division.size(); i++) {
            for (auto node : division[i]) {
                componentsOfNode[node].push_back(i);
            }
        }

        for (auto& components : componentsOfNode) {
            assert(components.size() > 0);
            assert(components.size() < 4);
        }

        for (const auto& [u, v] : Graph.edgeRange()) {
            bool hasCommon = std::any_of(componentsOfNode[u].begin(), componentsOfNode[u].end(), [&](int x) {
                return std::find(componentsOfNode[v].begin(), componentsOfNode[v].end(), x) != componentsOfNode[v].end();
                });
            assert(hasCommon && "edge must be in at least one region fully");
        }
    }

    void printGraph(const NetworKit::Graph& graph) {
        std::cout << "Graph structure" << std::endl;
        for (auto edge : graph.edgeWeightRange()) {
            std::cout << edge.u << " - " << edge.v << " weight: " << edge.weight << std::endl;
        }
    }

    void print_division(nodeSubsets_t division) {
        std::cout << "Print Division" << std::endl;
        for (auto& div : division) {
            for (auto n : div) {
                std::cout << n << " ";
            }std::cout << std::endl;
        }
    }
}
