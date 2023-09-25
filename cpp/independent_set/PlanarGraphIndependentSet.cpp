/*
 * PlanarGraphIndependentSet.cpp
 *
 *  Created on: 06.04.2023
 *      Author: Mikołaj Twaróg
 *      Ported by: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include <independent_set/PlanarGraphIndependentSet.hpp>

#include <map>
#include <optional>
#include <tuple>

#include <boost/graph/adjacency_list.hpp>

#include <techniques/baker/bakers_technique.hpp>

#include <techniques/baker/BakerKOuterplanar.hpp>
#include <techniques/baker/Bodlaender.hpp>
#include <techniques/baker/ProblemNode.hpp>

auto get_embedding(Graph& g) {
    PlanarEmbedding embedding(boost::num_vertices(g));
    boost::property_map<Graph, boost::edge_index_t>::type e_index = get(boost::edge_index, g);
    boost::graph_traits<Graph>::edges_size_type edge_count = 0;
    boost::graph_traits<Graph>::edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
        put(e_index, *ei, edge_count++);
    }
    boost::boyer_myrvold_planarity_test(g, &embedding[0]);
    std::map<std::pair<int, int>, std::vector<int>> faces;
    std::vector<std::vector<int>> vertices_in_face;
    face_getter visitor(faces, vertices_in_face);
    level_face_traversal2(embedding, visitor);
    return std::make_pair(embedding, vertices_in_face[0]);
}

namespace Koala {

BakerPlanarGraphIndependentSet::BakerPlanarGraphIndependentSet(
    NetworKit::Graph &graph, double epsilon)
        : PlanarGraphIndependentSet(graph), epsilon(epsilon) { }

void BakerPlanarGraphIndependentSet::run() {
    Graph G(graph->numberOfNodes());
    graph->forEdges([&](NetworKit::node u, NetworKit::node v) {
        boost::add_edge(u, v, G);
    });
    auto [embedding, outer_face] = get_embedding(G);
    int result = baker<IndependentSetNode>(G, *graph, embedding, outer_face);
    for (NetworKit::node u = 0; u < result; u++) {
        independentSet.insert(u);
    }
    hasRun = true;
}

BodlaenderPlanarGraphIndependentSet::BodlaenderPlanarGraphIndependentSet(
    NetworKit::Graph &graph, double epsilon)
        : PlanarGraphIndependentSet(graph), epsilon(epsilon) { }

void BodlaenderPlanarGraphIndependentSet::run() {
    Graph G(graph->numberOfNodes());
    graph->forEdges([&](NetworKit::node u, NetworKit::node v) {
        boost::add_edge(u, v, G);
    });
    auto [embedding, outer_face] = get_embedding(G);
    int result = bodlaender<IndependentSetNode>(G, *graph, embedding, outer_face);
    for (NetworKit::node u = 0; u < result; u++) {
        independentSet.insert(u);
    }
    hasRun = true;
}

}  /* namespace Koala */
