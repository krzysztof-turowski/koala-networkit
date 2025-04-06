#include <planar_sssp/PlanarUtils.hpp> 

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>
#include <boost/graph/biconnected_components.hpp>
#include <boost/graph/make_biconnected_planar.hpp>
#include <boost/graph/make_maximal_planar.hpp>
#include <set>
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


    NetworKit::Graph convertBoostToNetworKit(BoostGraph G, std::vector<NetworKit::node> nodeMap) {
        NetworKit::Graph result(num_vertices(G));

        for (auto [ei, ei_end] = edges(G); ei != ei_end; ++ei) {
            auto u = source(*ei, G);
            auto v = target(*ei, G);

            result.addEdge(nodeMap[u], nodeMap[v]);
        }
        return result;
    }

    void printEmbeding(const embedding_storage_t& embedding_storage,
        const BoostGraph& boostGraph) {
        for (size_t i = 0; i < embedding_storage.size(); ++i) {
            std::cout << "Vertex " << i << ": ";
            for (const auto& edge : embedding_storage[i]) {
                if (source(edge, boostGraph) == i)
                    std::cout << target(edge, boostGraph) << " ";
                else
                    std::cout << source(edge, boostGraph) << " ";
            }
            std::cout << "\n";
        }
    }

    planar_embedding_t findPlanarEmbeding(const NetworKit::Graph& G) {
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
                result[nodeMap[i]].push_back(nodeMap[node]);
            }
        }

        printEmbeding(embedding_storage, boostGraph);
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

        return convertBoostToNetworKit(g, nodeMap);
    }

}
