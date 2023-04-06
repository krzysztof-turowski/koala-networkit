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

#include "baker/bakers_technique.hpp"
#include "baker/level_face_traversal.hpp"
#include "baker/baker-k-outer-planar.hpp"
#include "baker/visitors.hpp"
#include "baker/problems.hpp"
#include "baker/name_levels.hpp"
#include "bodlaender/bodlaender_impl.hpp"

int independent_set_(Graph& g) {
    const int n = num_vertices(g);
    boost::dynamic_bitset<> s(n, 0);
    int mx = 0, last = n;
    while (s.count() != n) {
        graph_traits<Graph>::edge_iterator ei, ei_end;
        bool res = true;
        for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
            Edge e = *ei;
            if (s[e.m_source] && s[e.m_target]) {
                res = false;
                break;
            }
        }
        if (res) {
            int ones = s.count();
            mx = std::max(mx, ones);
        }
        for(int i = s.size() - 1; i >= 0; --i) {
            if ((s[i] ^= 0x1) == 0x1) {
                break;
            }
        }
    }
    return mx;
}

void get_embedding(Graph& g, PlanarEmbedding& embedding, std::vector<int>& outer_face) {
    boost::property_map<Graph, edge_index_t>::type e_index = get(edge_index, g);
    boost::graph_traits<Graph>::edges_size_type edge_count = 0;
    boost::graph_traits<Graph>::edge_iterator ei, ei_end;
    for(boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
        put(e_index, *ei, edge_count++);
    boost::boyer_myrvold_planarity_test(g, &embedding[0]);
    std::map<boost::graph_traits<Graph>::edge_descriptor, std::vector<int> > faces;
    std::vector<std::vector<int> > vertices_in_face;
    face_getter<Edge> my_vis(&faces, vertices_in_face);
    level_face_traversal<Graph>(embedding, my_vis);
    for (int v : vertices_in_face[0]) {
        outer_face.push_back(v);
    }
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

    PlanarEmbedding embedding(boost::num_vertices(G));
    std::vector<int> outer_face;
    get_embedding(G, embedding, outer_face);
    int result = bakers_technique(G, embedding, outer_face, 3, Baker, is);
    // int result = baker<independent_set>(G, embedding, outer_face);
    assert(result == independent_set_(G));
    hasRun = true;
}

}  /* namespace Koala */
