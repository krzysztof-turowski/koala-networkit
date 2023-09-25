/*
 * Visitors.hpp
 *
 *  Created on: 08.02.2021
 *      Author: Mikołaj Twaróg
 *      Ported by: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/planar_face_traversal.hpp>
#include <boost/graph/subgraph.hpp>

#include <structures/CyclicVector.hpp>
#include <techniques/baker/PlanarTree.hpp>

typedef boost::subgraph<boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::property<boost::vertex_index_t, int>, boost::property<boost::edge_index_t, int>>> Graph;
typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
typedef boost::graph_traits<Graph>::edge_descriptor Edge;
typedef std::vector<cyclic_vector<Edge>> PlanarEmbedding;
typedef std::vector<cyclic_vector<NetworKit::Edge>> PlanarEmbedding2;

struct face_getter : public boost::planar_face_traversal_visitor {
    std::map<std::pair<int, int>, std::vector<int>> &faces;
    std::vector<std::vector<int>> &vertices_in_face;
    int current_face = 0;

    face_getter(std::map<std::pair<int, int>, std::vector<int>> &f, std::vector<std::vector<int>> &o) : faces(f), vertices_in_face(o) { }

    void begin_face() {
        vertices_in_face.emplace_back();
    }

    void end_face() {
        current_face++;
    }

    void next_edge(Edge e) {
        faces[std::minmax(e.m_source, e.m_target)].push_back(current_face);
    }

    void next_edge(NetworKit::Edge e) {
        faces[std::minmax(e.u, e.v)].push_back(current_face);
    }

    void next_vertex(NetworKit::node v) {
        vertices_in_face[current_face].push_back(v);
    }
};

template <typename Problem>
struct tree_builder : public boost::planar_face_traversal_visitor {
    const std::map<std::pair<int, int>, std::vector<int>> &faces;
    PlanarTree<Problem> &tree;
    const Graph &graph;
    int current_face = 0;
    int last_vertex;

    tree_builder(std::map<std::pair<int, int>, std::vector<int>> &f, PlanarTree<Problem> &t, Graph &g) : faces(f), tree(t), graph(g) { }

    void end_face() {
        current_face++;
    }

    void next_vertex(NetworKit::node v) {
        tree[current_face].face.push_back(graph.local_to_global(v));
        last_vertex = graph.local_to_global(v);
    }

    void next_edge(Edge e) {
        if(current_face != tree.outer_face) {
            const auto &f = faces.find(std::minmax(e.m_source, e.m_target))->second;
            int neighbor = f[0] == current_face ? f[1] : f[0];
            if (neighbor == tree.outer_face) {
                int last = tree.size();
                tree.emplace_back();
                tree[current_face].children.push_back(last);
                tree[last].parent = current_face;
                tree[last].label.second = last_vertex;
                tree[last].label.first = last_vertex == e.m_source ? e.m_target : e.m_source;
            } else {
                tree[current_face].children.push_back(neighbor);
            }
        }
    }

    void next_edge(NetworKit::Edge e) {
        if(current_face != tree.outer_face) {
            const auto &f = faces.find(std::minmax(e.u, e.v))->second;
            int neighbor = f[0] == current_face ? f[1] : f[0];
            if (neighbor == tree.outer_face) {
                int last = tree.size();
                tree.emplace_back();
                tree[current_face].children.push_back(last);
                tree[last].parent = current_face;
                tree[last].label.second = last_vertex;
                tree[last].label.first = last_vertex == e.u ? e.v : e.u;
            } else {
                tree[current_face].children.push_back(neighbor);
            }
        }
    }
};
