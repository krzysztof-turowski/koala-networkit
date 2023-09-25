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

template <typename Visitor>
inline void level_face_traversal2(PlanarEmbedding &embedding, Visitor &visitor) {
    std::unordered_map<NetworKit::Edge, std::unordered_map<NetworKit::node, NetworKit::Edge>> next_edge;
    std::unordered_map<NetworKit::Edge, std::unordered_set<NetworKit::node>> visited;
    std::vector<NetworKit::Edge> edges_cache;

    visitor.begin_traversal();
    for (int v = 0; v < embedding.size(); v++) {
        auto& edges = embedding[v];
        for (int i = 0; i < edges.size(); i++) {
            auto e = NetworKit::Edge(edges[i].m_source, edges[i].m_target);
            edges_cache.push_back(e);
            next_edge[e][v] = NetworKit::Edge(edges[i + 1].m_source, edges[i + 1].m_target);
        }
    }
    for (auto e : edges_cache) {
        for (auto v : {e.u, e.v}) {
            auto e_visited_found = visited[e].find(v);
            if (e_visited_found == visited[e].end()) {
                visitor.begin_face();
            }
            while (visited[e].find(v) == visited[e].end()) {
                visitor.next_vertex(v);
                visitor.next_edge(e);
                visited[e].insert(v);
                v = e.u == v ? e.v : e.u;
                e = next_edge[e][v];
            }
            if (e_visited_found == visited[e].end()) {
                visitor.end_face();
            }
        }
    }
}

template <typename Visitor>
inline void level_face_traversal(PlanarEmbedding2 &embedding, Visitor &visitor) {
    std::unordered_map<NetworKit::Edge, std::unordered_map<NetworKit::node, NetworKit::Edge>> next_edge;
    std::unordered_map<NetworKit::Edge, std::unordered_set<NetworKit::node>> visited;
    std::vector<NetworKit::Edge> edges_cache;

    visitor.begin_traversal();
    for (int v = 0; v < embedding.size(); v++) {
        auto& edges = embedding[v];
        for (int i = 0; i < edges.size(); i++) {
            auto &e = edges[i];
            edges_cache.push_back(e);
            next_edge[e][v] = edges[i + 1];
        }
    }
    for (auto e : edges_cache) {
        for (auto v : {e.u, e.v}) {
            auto e_visited_found = visited[e].find(v);
            if (e_visited_found == visited[e].end()) {
                visitor.begin_face();
            }
            while (visited[e].find(v) == visited[e].end()) {
                visitor.next_vertex(v);
                visitor.next_edge(e);
                visited[e].insert(v);
                v = e.u == v ? e.v : e.u;
                e = next_edge[e][v];
            }
            if (e_visited_found == visited[e].end()) {
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

int get_edge_it(int v, int w, const PlanarEmbedding& embedding) {
    auto element = std::find_if(
        embedding[v].begin(), embedding[v].end(), [&](auto &e) {
            return (e.m_source == v && e.m_target == w) || (e.m_source == w && e.m_target == v);
        });
    return element != embedding[v].end() ? std::distance(embedding[v].begin(), element) : -1;
}

int get_edge_it(Edge e, int v, const PlanarEmbedding& embedding) {
    int w = e.m_source == v ? e.m_target : e.m_source;
    return get_edge_it(v, w, embedding);
}

int get_edge_it(NetworKit::Edge e, int v, const PlanarEmbedding& embedding) {
    int w = e.u == v ? e.v : e.u;
    return get_edge_it(v, w, embedding);
}

int get_edge_it(NetworKit::node v, NetworKit::node w, const PlanarEmbedding2 &embedding) {
    auto element = std::find_if(
        embedding[v].begin(), embedding[v].end(), [&](auto &e) {
            return (e.u == v && e.v == w) || (e.u == w && e.v == v);
        });
    return element != embedding[v].end() ? std::distance(embedding[v].begin(), element) : -1;
}

int get_edge_it(NetworKit::Edge e, NetworKit::node v, const PlanarEmbedding2 &embedding) {
    int w = e.u == v ? e.v : e.u;
    return get_edge_it(v, w, embedding);
}
