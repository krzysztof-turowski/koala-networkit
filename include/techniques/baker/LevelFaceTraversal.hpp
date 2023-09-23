/*
 * LevelFaceTraversal.hpp
 *
 *  Created on: 04.02.2021
 *      Author: Mikołaj Twaróg
 *      Ported by: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

#include <vector>
#include <set>
#include <map>

#include <boost/graph/boyer_myrvold_planar_test.hpp>

#include <techniques/baker/Visitors.hpp>

template <typename Graph, typename Visitor, typename Container>
inline void level_face_traversal(Container& embedding, Visitor& visitor) {
    std::map<Edge, std::map<Vertex, Edge>> next_edge;
    std::map<Edge, std::set<Vertex>> visited;

    visitor.begin_traversal();
    std::vector<Edge> edges_cache;

    for (int v = 0; v < embedding.size(); v++) {
        auto& edges = embedding[v];
        for (int i = 0; i < edges.size(); i++) {
            Edge e = edges[i];
            edges_cache.push_back(e);
            next_edge[e][v] = edges[(i + 1) % edges.size()];
        }
    }

    for (auto e : edges_cache) {
        for (auto v : {e.m_source, e.m_target}) {
            std::set<Vertex>& e_visited = visited[e];
            auto e_visited_found = e_visited.find(v);
            if (e_visited_found == visited[e].end()) {
                visitor.begin_face();
            }
            while (visited[e].find(v) == visited[e].end()) {
                visitor.next_vertex(v);
                visitor.next_edge(e);
                visited[e].insert(v);
                v = e.m_source == v ? e.m_target : e.m_source;
                e = next_edge[e][v];
            }
            if (e_visited_found == e_visited.end()) {
                visitor.end_face();
            }
        }
    }
}

template<typename Graph>
bool check_for_edge(int x, int y, Graph& g, std::set<std::pair<int, int>>& added_edges) {
    if (x == y) {
        return false;
    }
    bool res = false;
    typename boost::graph_traits<Graph>::out_edge_iterator ei, ei_end;
    for(boost::tie(ei, ei_end) = out_edges(x, g); ei != ei_end; ++ei) {
        if (ei->m_source == y || ei->m_target == y) {
            res = true;
            break;
        }
    }
    std::pair<int, int> e(x, y);
    if (!res || added_edges.find(e) != added_edges.end()) {
        return false;
    }
    std::swap(e.first, e.second);
    return added_edges.find(e) == added_edges.end();
}

int get_edge_it(int v, int w, PlanarEmbedding& embedding) {
    auto element = std::find_if(
        embedding[v].begin(), embedding[v].end(), [&](auto &e) {
            return (e.m_source == v && e.m_target == w) || (e.m_source == w && e.m_target == v);
        });
    return element != embedding[v].end() ? std::distance(embedding[v].begin(), element) : -1;
}

int get_edge_it(Edge e, int v, PlanarEmbedding& embedding) {
    int w = e.m_source == v ? e.m_target : e.m_source;
    return get_edge_it(v, w, embedding);
}

auto get_embedding(const NetworKit::Graph &graph) {
    Graph g(graph.numberOfNodes());
    graph.forEdges([&](NetworKit::node u, NetworKit::node v) {
        boost::add_edge(u, v, g);
    });
    std::vector<cyclic_vector<Edge>> embedding(boost::num_vertices(g));
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
    level_face_traversal<Graph>(embedding, visitor);
    /*PlanarEmbedding out_embedding(boost::num_vertices(g));
    for (int i = 0; i < embedding.size(); i++) {
        out_embedding.reserve(embedding[i].size());
        for (auto e : embedding[i]) {
            out_embedding.emplace(e.m_source, e.m_target);
        }
    }*/
    return std::make_tuple(g, embedding, vertices_in_face[0]);
}
