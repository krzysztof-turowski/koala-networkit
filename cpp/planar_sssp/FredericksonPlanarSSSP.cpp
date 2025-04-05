#include "planar_sssp/FredericksonPlanarSSSP.hpp"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>
#include <set>
#include <unordered_map>
#include <utility>

using namespace boost;

// Boost graph typedefs
using BoostGraph =
    adjacency_list<vecS, vecS, undirectedS, property<vertex_index_t, int>>;

// Boost planar embedding typedefs based on
// https://www.boost.org/doc/libs/1_53_0/libs/graph/example/straight_line_drawing.cpp
using embedding_storage_t =
    std::vector<std::vector<graph_traits<BoostGraph>::edge_descriptor>>;
using embedding_t = boost::iterator_property_map<
    embedding_storage_t::iterator,
    property_map<BoostGraph, vertex_index_t>::type>;

namespace Koala
{

    BoostGraph convertNetworKitToBoost(const NetworKit::Graph &netGraph)
    {
        BoostGraph boostGraph(netGraph.numberOfNodes());
        netGraph.forEdges([&](NetworKit::node u, NetworKit::node v)
                          { add_edge(u, v, boostGraph); });
        return boostGraph;
    }

    void printEmbeding(const embedding_storage_t &embedding_storage,
                       const BoostGraph &boostGraph)
    {
        for (size_t i = 0; i < embedding_storage.size(); ++i)
        {
            std::cout << "Vertex " << i << ": ";
            for (const auto &edge : embedding_storage[i])
            {
                if (source(edge, boostGraph) == i)
                    std::cout << target(edge, boostGraph) << " ";
                else
                    std::cout << source(edge, boostGraph) << " ";
            }
            std::cout << "\n";
        }
    }

    void printGraph(const NetworKit::Graph &graph)
    {
        for (NetworKit::node u = 0; u < graph.numberOfNodes(); ++u)
        {
            for (NetworKit::node v : graph.neighborRange(u))
            {
                if (u < v)
                {
                    std::cout << u << " - " << v << std::endl;
                }
            }
        }
    }

    NetworKit::Graph convertToMaxDegree3(NetworKit::Graph &G)
    {
        BoostGraph boostGraph = convertNetworKitToBoost(G);
        embedding_storage_t embedding_storage(num_vertices(boostGraph));
        embedding_t embedding(embedding_storage.begin(),
                              get(vertex_index, boostGraph));

        if (!boyer_myrvold_planarity_test(
                boyer_myrvold_params::graph = boostGraph,
                boyer_myrvold_params::embedding = embedding))
        {
            throw std::runtime_error("Graph have to be planar!");
            return;
        }

        printEmbeding(embedding_storage, boostGraph);

        NetworKit::count result_size = G.numberOfNodes();
        for (size_t i = 0; i < embedding_storage.size(); ++i)
        {
            if (embedding_storage[i].size() > 3)
            {
                result_size += embedding_storage[i].size() - 1;
            }
        }

        if (result_size == G.numberOfNodes())
        {
            return G;
        }

        std::vector<std::unordered_map<size_t, size_t>> vertexMaps(
            G.numberOfNodes(), std::unordered_map<size_t, size_t>());
        std::vector<std::vector<size_t>> embeddingV(G.numberOfNodes(),
                                                    std::vector<size_t>());

        for (size_t i = 0; i < embedding_storage.size(); ++i)
        {
            for (auto &edge : embedding_storage[i])
            {
                size_t t = source(edge, boostGraph) == i ? target(edge, boostGraph)
                                                         : source(edge, boostGraph);
                embeddingV[i].push_back(t);
            }
        }

        size_t emptyVertex = G.numberOfNodes();
        NetworKit::Graph result(result_size);

        for (size_t i = 0; i < embeddingV.size(); ++i)
        {
            if (embeddingV[i].size() > 3)
            {
                size_t t = embeddingV[i][0];
                vertexMaps[i][t] = i;
                vertexMaps[i][emptyVertex] = 0;

                result.addEdge(i, emptyVertex, 0);
                for (int j = 1; j < embedding_storage[i].size(); ++j)
                {
                    t = embeddingV[i][j];
                    vertexMaps[i][t] = emptyVertex;
                    if (j < embedding_storage[i].size() - 1)
                    {
                        result.addEdge(emptyVertex, emptyVertex + 1, 0);
                    }
                    emptyVertex++;
                }

                result.addEdge(emptyVertex - 1, i, 0);
            }
        }

        G.forEdges(
            [&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w)
            {
                size_t s = v;
                size_t t = u;
                if (vertexMaps[u].find(v) != vertexMaps[u].end())
                {
                    t = vertexMaps[u][v];
                }
                if (vertexMaps[v].find(u) != vertexMaps[v].end())
                {
                    s = vertexMaps[v][u];
                }
                result.addEdge(s, t, w);
            });

        return result;
    }

    void FredericksonPlanarSSSP::run()
    {
        normal_graph = convertToMaxDegree3(graph);

        printGraph(normal_graph);
        hasRun = true;
        distanceToTarget = 5;
        return;
    }

} /* namespace Koala */
