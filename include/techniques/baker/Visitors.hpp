/*
 * Visitors.hpp
 *
 *  Created on: 08.02.2021
 *      Author: Mikołaj Twaróg
 *      Ported by: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

template <typename Edge>
struct face_getter : public boost::planar_face_traversal_visitor {
    std::map<std::pair<int, int>, std::vector<int>> &faces;
    std::vector<std::vector<int>> &vertices_in_face;
    int current_face = 0;

    face_getter(std::map<std::pair<int, int>, std::vector<int>> &f, std::vector<std::vector<int>> &o)
        : faces(f), vertices_in_face(o) { }

    void begin_face() {
        vertices_in_face.emplace_back();
    }

    void end_face() {
        current_face++;
    }

    void next_edge(Edge e) {
        faces[std::minmax(e.m_source, e.m_target)].push_back(current_face);
    }

    template <typename Vertex>
    void next_vertex(Vertex v) {
        vertices_in_face[current_face].push_back(v);
    }
};

template <typename Edge, typename Problem, typename PlanarEmbedding>
struct tree_builder : public boost::planar_face_traversal_visitor {
    std::map<std::pair<int, int>, std::vector<int>> &faces;
    PlanarTree<Problem> &tree;
    int current_face = 0;
    int last_vertex;
    Graph graph;

    tree_builder(std::map<std::pair<int, int>, std::vector<int>> &f, PlanarTree<Problem> &t, Graph &g)
            : faces(f), tree(t), graph(g) { }

    void end_face() {
        current_face++;
    }

    template <typename Vertex>
    void next_vertex(Vertex v) {
        tree[current_face].face.push_back(graph.local_to_global(v));
        last_vertex = graph.local_to_global(v);
    }

    void next_edge(Edge e) {
        if(current_face != tree.outer_face) {
            auto key = std::minmax(e.m_source, e.m_target);
            int neighbor = faces[key][0] == current_face ? faces[key][1] : faces[key][0];
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
};
