//
// Created by mikolajtwarog on 2021-02-04.
//

#ifndef TECHNIKA_BAKER_LEVEL_FACE_TRAVERSAL_HPP
#define TECHNIKA_BAKER_LEVEL_FACE_TRAVERSAL_HPP

#include <vector>
#include <set>
#include <map>
#include <boost/next_prior.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

template < typename Graph, typename Visitor , typename Container>
inline void level_face_traversal(Container& embedding, Visitor& visitor) {
    typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_t;
    typedef typename boost::graph_traits<Graph>::edge_descriptor edge_t;

    using h = std::hash<int>;
    auto hash = [](const edge_t& n){return ((17 * 31 + h()(n.m_source)) * 31 + h()(n.m_target)) * 31;};
    auto equal = [](const edge_t & l, const edge_t& r){return l.m_source == r.m_source && l.m_target == r.m_target && l.m_eproperty == r.m_eproperty;};
    std::map<edge_t, std::map<vertex_t, edge_t> > next_edge;
    std::map<edge_t, std::set<vertex_t> > visited;

    visitor.begin_traversal();
    std::vector<edge_t> edges_cache;

    for (int v = 0; v < embedding.size(); v++) {
        auto& edges = embedding[v];

        for (int i = 0; i < edges.size(); i++) {
            edge_t e = edges[i];

            edges_cache.push_back(e);

            next_edge[e][v] = edges[(i+1)%edges.size()];
        }
    }

    std::vector<vertex_t> vertices_in_edge;

    for (auto e : edges_cache) {
        vertices_in_edge.clear();
        vertices_in_edge.push_back(e.m_source);
        vertices_in_edge.push_back(e.m_target);

        for (auto v : vertices_in_edge) {
            std::set<vertex_t>& e_visited = visited[e];
            typename std::set<vertex_t>::iterator e_visited_found
                = e_visited.find(v);

            if (e_visited_found == visited[e].end())
                visitor.begin_face();

            while (visited[e].find(v) == visited[e].end())
            {
                visitor.next_vertex(v);
                visitor.next_edge(e);
                visited[e].insert(v);
                v = e.m_source == v ? e.m_target : e.m_source;
                e = next_edge[e][v];
            }

            if (e_visited_found == e_visited.end())
                visitor.end_face();
        }
    }
}

template<typename Graph>
bool check_for_edge(int x, int y, Graph& g, std::set<std::pair<int, int> >& added_edges) {
    if (x == y) {
        return false;
    }

    std::pair<int, int> e(x, y);

    bool res = false;
    typename graph_traits<Graph>::out_edge_iterator ei, ei_end;
    for(boost::tie(ei, ei_end) = out_edges(x, g); ei != ei_end; ++ei) {
        if (ei->m_source == y || ei->m_target == y) {
            res = true;
            break;
        }
    }

    if (!res || added_edges.find(e) != added_edges.end()) {
        return false;
    }

    std::swap(e.first, e.second);

    return added_edges.find(e) == added_edges.end();
}

int get_edge_it(int v, int w, PlanarEmbedding& embedding) {
    for (int i = 0; i < embedding[v].size(); i++) {
        if ((embedding[v][i].m_source == v
             && embedding[v][i].m_target == w)
            || (embedding[v][i].m_source == w
                && embedding[v][i].m_target == v)) {
            return i;
        }
    }

    return -1;
}

int get_edge_it(Edge e, int v, PlanarEmbedding& embedding) {
    int w = e.m_source == v ? e.m_target : e.m_source;
    return get_edge_it(v, w, embedding);
}

#endif //TECHNIKA_BAKER_LEVEL_FACE_TRAVERSAL_HPP
