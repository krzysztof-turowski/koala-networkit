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

#include <boost/graph/graph_traits.hpp>
#include <boost/ref.hpp>

typedef boost::subgraph<boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::property<boost::vertex_index_t, int>, boost::property<boost::edge_index_t, int>>> Graph;
typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
typedef boost::graph_traits<Graph>::edge_descriptor Edge;
typedef std::vector<cyclic_vector<Edge>> PlanarEmbedding;

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
            next_edge[e][v] = edges[(i + 1)%edges.size()];
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
