//
// Created by mikolajtwarog on 2021-02-08.
//

#ifndef TECHNIKA_BAKER_BAKER_K_OUTER_PLANAR_HPP

#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/ref.hpp>
#include <vector>
#include <climits>
#include <queue>

#include <boost/graph/planar_canonical_ordering.hpp>
#include <boost/graph/is_straight_line_drawing.hpp>
#include <boost/graph/chrobak_payne_drawing.hpp>
#include <boost/graph/planar_face_traversal.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>
#include <boost/graph/make_maximal_planar.hpp>
#include <boost/graph/subgraph.hpp>
#include "../utils/cyclic_vector.hpp"

using namespace boost;

typedef subgraph<adjacency_list
        <
                vecS,
                vecS,
                undirectedS,
                property<vertex_index_t, int>,
                property<edge_index_t, int>
        > >
        Graph;

typedef std::vector<cyclic_vector< graph_traits<Graph>::edge_descriptor > > PlanarEmbedding;

typedef graph_traits<Graph>::edge_descriptor Edge;

typedef graph_traits<Graph>::vertex_descriptor Vertex;

#include "../utils/level_face_traversal.hpp"
#include "problems.hpp"
#include "../utils/visitors.hpp"
#include "../utils/name_levels.hpp"

template<typename Problem>
class baker_impl {
    Graph g;
    PlanarEmbedding embedding;
    std::vector<int> outer_face;
    ::tree<Problem> tree;
    std::vector<int> vertex_level;
    std::set< std::pair<int, int> > added_edges;

    void triangulate(std::vector<int>& face, std::vector<int>& component, int turn) {
        if (component.size() == 1) {
            int c = component[0];
            embedding[c].clear();

            for (int i = 0; i < face.size(); i++) {
                int v = face[i];

                embedding[c].emplace_back(v, c, &c);

                if (!edge(v, c, g).second) {
                    Edge new_e = add_edge(c, v, g).first;
                    int next = face[(i + 1) % face.size()];
                    embedding[v].insert(embedding[v].begin() + get_edge_it(v, next, embedding) + 1,
                                        new_e);
                    added_edges.emplace(v, c);
                }
            }
            return;
        }

        int level = vertex_level[face[0]];
        Edge connecting_e;
        int starting_v;
        int starting_v_it = -1;
        int c_v;

        for (int i = 0; i < face.size(); i++) {
            int v = face[i];
            for (Edge e : embedding[v]) {
                int neighbour = e.m_source == v ? e.m_target : e.m_source;
                if (std::find(component.begin(), component.end(), neighbour) != component.end()) {
                    connecting_e = e;
                    starting_v = v;
                    starting_v_it = i;
                    c_v = neighbour;
                    break;
                }
            }
            if (starting_v_it > -1)
                break;
        }

        std::rotate(face.begin(), face.begin() + starting_v_it, face.end());
        std::rotate(component.begin(), std::find(component.begin(), component.end(), c_v), component.end());

        int face_it = 0;
        int comp_it = 0;
        int curr_e_it_face = get_edge_it(connecting_e, starting_v, embedding);
        int curr_e_it_comp = get_edge_it(connecting_e, c_v, embedding);

        std::set<int> comp_vis;
        std::set<int> face_vis;

        int comp = 0;

        while (comp < component.size() || face_vis.size() < face.size()) {
            bool res = false;
            int comp_curr = component[comp_it];
            int face_curr = face[face_it];
            int last = -1;

            face_vis.insert(face_curr);

            int first_e = get_edge_it(comp_curr, component[(comp_it - 1 + component.size()) % component.size()], embedding);

            int second_e = get_edge_it(comp_curr, component[(comp_it + 1) % component.size()], embedding);

            for (int i = (first_e + 1) % embedding[comp_curr].size(); i != second_e;
                 i = (i + 1) % embedding[comp_curr].size()) {
                 if (i == curr_e_it_comp) {
                     res = true;
                     break;
                 }
            }

            if (!res) {
                curr_e_it_comp = first_e;
            }

            res = false;

            for (int i = (curr_e_it_comp + 1) % embedding[comp_curr].size(); i != second_e;
            i = (i + 1) % embedding[comp_curr].size()) {
                Edge& e = embedding[comp_curr][i];
                int neighbour = e.m_source == comp_curr ? e.m_target : e.m_source;
                if (vertex_level[neighbour] == level && neighbour != face_curr) {
                    res = true;
                    last = neighbour;
                    break;
                }
            }

            if (res)  {
                for (int i = (face_it + 1) % face.size(); ; i = (i + 1) % face.size()) {
                    if (face[i] == last) {
                        face_it = i;
                        curr_e_it_face = (get_edge_it(last, face[(i - 1 + face.size()) % face.size()], embedding)
                                - 1 + embedding[last].size()) % embedding[last].size();
                        curr_e_it_comp = (curr_e_it_comp + 1) % embedding[comp_curr].size();
                        break;
                    }

                    face_vis.insert(face[i]);
                    if (!edge(comp_curr, face[i], g).second) {
                        added_edges.emplace(face[i], comp_curr);
                    }

                    curr_e_it_comp = (curr_e_it_comp + 1) % embedding[comp_curr].size();
                    Edge e(face[i - 1], face[i], nullptr);

                    int e_it;
                    if (i == face_it) {
                        e_it = curr_e_it_face;
                    } else {
                        e_it = get_edge_it(e, face[i], embedding);
                    }


                    if (embedding[face[i]][e_it].m_source != comp_curr
                        && embedding[face[i]][e_it].m_target != comp_curr) {
                        Edge new_e = add_edge(face[i], comp_curr, g).first;

                        embedding[face[i]].insert(embedding[face[i]].begin() + e_it, new_e);
                        embedding[comp_curr].insert(embedding[comp_curr].begin() +
                                                    curr_e_it_comp, new_e);
                    }
                }
            } else {
                comp_vis.insert(comp_curr);
                comp++;
                comp_it = (comp_it + 1) % component.size();
                curr_e_it_comp = get_edge_it(component[comp_it], comp_curr, embedding);
                comp_curr = component[comp_it];

                if (!edge(comp_curr, face_curr, g).second) {
                    added_edges.emplace(face_curr, comp_curr);
                }

                int temp_e_it_face = (curr_e_it_face - 1 + embedding[face_curr].size()) % embedding[face_curr].size();
                int temp_e_it_comp = (curr_e_it_comp + 1) % embedding[comp_curr].size();

                if (embedding[comp_curr][temp_e_it_comp].m_source != face_curr
                && embedding[comp_curr][temp_e_it_comp].m_target != face_curr) {

                    Edge new_e = add_edge(face_curr, comp_curr, g).first;

                    embedding[face_curr].insert(embedding[face_curr].begin() + curr_e_it_face,
                                                new_e);
                    embedding[comp_curr].insert(embedding[comp_curr].begin() +
                                                (curr_e_it_comp + 1),
                                                new_e);
                } else {
                    Edge& e_con = embedding[comp_curr][temp_e_it_comp];
                    curr_e_it_face = (curr_e_it_face - 1 + embedding[face_curr].size()) % embedding[face_curr].size();
                    for (; curr_e_it_face < embedding[face_curr].size(); curr_e_it_face++) {
                        if (embedding[face_curr][curr_e_it_face] == e_con) {
                            break;
                        }
                    }
                }
                curr_e_it_comp = (curr_e_it_comp + 1) % embedding[comp_curr].size();
            }
        }
    }

    int find_third(int one, int two, std::vector<int>& component) {
        auto& one_edges = embedding[one];
        auto& two_edges = embedding[two];

        int level = vertex_level[one];

        if (one == two) {
            for (auto e_i : one_edges) {
                int target_i = e_i.m_source == one ? e_i.m_target : e_i.m_source;

                if (vertex_level[target_i] == level + 1) {
                    return target_i;
                }
            }
        }

        int edge_it = get_edge_it(one, two, embedding);
        int cos = edge_it - one_edges.size();
        int third;
        int connecting_e_it = -1;
        for (int i = edge_it - 1; i > cos; i--) {
            third = one_edges[i].m_source == one ? one_edges[i].m_target : one_edges[i].m_source;

            if (third != one) {
                Edge curr_e;
                for (int j = 0; j < embedding[third].size(); j++) {
                    curr_e = embedding[third][j];
                    if (curr_e == one_edges[i]) {
                        connecting_e_it = j;
                        break;
                    }
                }
                if (connecting_e_it == -1) {
                    connecting_e_it = get_edge_it(one_edges[i], third, embedding);
                }
                break;
            }
        }

        int starting_v = third;
        component.push_back(starting_v);

        Edge current_e;
        int current_v;
        bool res = false;
        level++;
        auto &edges = embedding[starting_v];
        for (int i = (connecting_e_it + 1) % edges.size(); i != connecting_e_it;
             i = (i + 1) % edges.size()) {
            Edge e = edges[i];
            int neighbour = e.m_source == starting_v ? e.m_target : e.m_source;

            if (vertex_level[neighbour] == level) {
                current_e = e;
                current_v = neighbour;
                res = true;
                break;
            }
        }
        if (!res) {
            return third;
        }

        while (true) {
            component.push_back(current_v);

            int temp_v = current_v;

            int current_e_it = get_edge_it(current_e, current_v, embedding);
            auto &edges2 = embedding[current_v];
            int i = (current_e_it + 1 + edges2.size()) % edges2.size();
            res = false;
            for (; i != current_e_it; i = (i + 1 + edges2.size()) % edges2.size()) {
                Edge e = edges2[i];
                int neighbour = e.m_source == current_v ? e.m_target : e.m_source;

                if (vertex_level[neighbour] == level) {
                    current_e = e;
                    current_v = neighbour;
                    res = true;
                    break;
                }
            }

            if (res == false) {
                Edge e = edges2[current_e_it];
                int neighbour = e.m_source == current_v ? e.m_target : e.m_source;

                if (vertex_level[neighbour] == level) {
                    current_e = e;
                    current_v = neighbour;
                }
            }

            if (temp_v == starting_v && current_v == component[1]) {
                break;
            }
        }

        component.pop_back();

        return third;
    }

    int find_dividing_points(Edge one, Edge two, std::set<int>& dividing_points) {
        int v = one.m_target;

        int level = vertex_level[v];

        int one_it = get_edge_it(one, v, embedding);
        int two_it = get_edge_it(two, v, embedding);

        auto& edges = embedding[v];

        for (int i = (one_it + 1) % edges.size(); i != two_it; i = (i + 1) % edges.size()) {
            Edge e = edges[i];
            if (vertex_level[e.m_source] == level - 1 || vertex_level[e.m_target] == level - 1) {
                 dividing_points.insert(vertex_level[e.m_source] == level - 1 ? e.m_source : e.m_target);
                 return vertex_level[e.m_source] == level - 1 ? e.m_source : e.m_target;
            }
        }

        return -1;
    }

    void get_leaves(::tree<Problem>& t, std::vector<int>& leaves, int v) {
        if (t[v].children.empty()) {
            leaves.push_back(v);
            return;
        }

        for (auto child : t[v].children) {
            get_leaves(t, leaves, child);
        }
    }

    void get_component(std::vector< std::vector<int> >& components, std::map<int, int>& vis,
                       std::vector< std::pair<int, int> >& v_in_c) {
        if (vertex_level[v_in_c[0].first] == 1) {
            for (int v : outer_face) {
                components[0].push_back(v);
            }

            return;
        }

        int comp_num = 0;

        for (int c = 0; c < v_in_c.size(); c++) {
            if (vis.find(v_in_c[c].first) != vis.end()) {
                continue;
            }
            components.emplace_back();

            int v = v_in_c[c].first;

            int level = vertex_level[v];
            int starting_v = v;
            int connecting_e_it = get_edge_it(v, v_in_c[c].second, embedding);

            vis[starting_v] = comp_num;
            components.back().push_back(starting_v);

            Edge current_e;
            int current_v;
            bool res = false;

            auto &edges = embedding[starting_v];
            for (int i = (connecting_e_it + 1) % edges.size(); i != connecting_e_it;
                 i = (i + 1) % edges.size()) {
                Edge e = edges[i];
                int neighbour = e.m_source == starting_v ? e.m_target : e.m_source;

                if (vertex_level[neighbour] == level) {
                    current_e = e;
                    current_v = neighbour;
                    res = true;
                    break;
                }
            }
            if (!res) {
                comp_num++;
                continue;
            }

            while (true) {
                vis[current_v] = comp_num;
                components.back().push_back(current_v);

                int temp_v = current_v;

                int current_e_it = get_edge_it(current_e, current_v, embedding);
                auto &edges2 = embedding[current_v];
                int i = (current_e_it + 1 + edges2.size()) % edges2.size();
                res = false;
                for (; i != current_e_it; i = (i + 1 + edges2.size()) % edges2.size()) {
                    Edge e = edges2[i];
                    int neighbour = e.m_source == current_v ? e.m_target : e.m_source;

                    if (vertex_level[neighbour] == level) {
                        current_e = e;
                        current_v = neighbour;
                        res = true;
                        break;
                    }
                }

                if (res == false) {
                    Edge e = edges2[current_e_it];
                    int neighbour = e.m_source == current_v ? e.m_target : e.m_source;

                    if (vertex_level[neighbour] == level) {
                        current_e = e;
                        current_v = neighbour;
                    }
                }

                if (temp_v == starting_v && current_v == components.back()[1]) {
                    break;
                }
            }

            components.back().pop_back();
            comp_num++;
        }
    }

    void make_connected(std::vector< std::vector<int> >& components, std::map<int, int>& vis,
                        std::vector< std::pair<int, int> >& v_in_c) {
        std::vector<bool> connected(components.size(), false);

        for (int i = 0; i < v_in_c.size(); i++) {
            int curr = v_in_c[i].first;
            int curr_con = v_in_c[i].second;
            int next = v_in_c[(i + 1) % v_in_c.size()].first;
            int next_con = v_in_c[(i + 1) % v_in_c.size()].second;

            connected[vis[curr]] = true;

            if (vis[curr] == vis[next] || edge(curr, next, g).second || connected[vis[next]]) {
                continue;
            }

            added_edges.insert(std::pair<int, int>(curr, next));
            Edge new_e = add_edge(curr, next, g).first;

            embedding[curr].insert(embedding[curr].begin() + get_edge_it(curr, curr_con, embedding) + 1, new_e);

            embedding[next].insert(embedding[next].begin() + get_edge_it(next, next_con, embedding), new_e);
        }
    }

    void check_for_components(const std::vector<int> &face, std::vector< std::pair<int, int> >& out) {
        if (face.size() == 2) {
            return;
        }

        int level = vertex_level[face[0]];

        int prev = face.size() - 1, curr = 0, next = 1;

        for (; curr < face.size(); prev = (prev + 1) % face.size(), curr++, next = (next + 1) % face.size()) {
            int one = get_edge_it(face[curr], face[prev], embedding);
            int two = get_edge_it(face[curr], face[next], embedding);
            for (int i = (one - 1 + embedding[face[curr]].size()) % embedding[face[curr]].size(); i != two;
            i = (i - 1 + embedding[face[curr]].size()) % embedding[face[curr]].size()) {
                Edge e = embedding[face[curr]][i];
                int target = e.m_source == face[curr] ? e.m_target : e.m_source;

                if (vertex_level[target] == level + 1) {
                    out.emplace_back(target, face[curr]);
                }
            }
        }

        return;
    }

    void root_tree_2(::tree<Problem> &t, int node) {
        if (t[node].children.empty())
            return;

        int parent_it = 0;
        auto child = t[node].children.begin();
        for (int i = 0; child != t[node].children.end(); child++, i++) {
            if (*child == t[node].parent) {
                t[node].children.erase(child--);
                parent_it = i;
            } else {
                t[*child].parent = node;
                root_tree_2(t, *child);
            }
        }

        int last = t[node].children.size() - 1;

        std::rotate(t[node].children.begin(),
                    t[node].children.begin() + parent_it,
                    t[node].children.end());

        t[node].label.first = t[t[node].children[0]].label.first;
        t[node].label.second = t[t[node].children[last]].label.second;
    }

    void root_tree_with_root(::tree<Problem> &t, int root) {
        int node;

        for (int i = 0; i < t.size(); i++) {
            for (int v : t[i].face) {
                if (v == root) {
                    node = i;
                    break;
                }
            }
        }

        t.root = node;

        root_tree_2(t, node);

        auto &children = t[node].children;

        for (int i = 0; i < children.size(); i++) {
            int child = children[i];

            if (t[child].label.first == root) {
                std::rotate(t[node].children.begin(),
                            t[node].children.begin() + i,
                            t[node].children.end());
                break;
            }
        }

        t[node].label.first = root;
        t[node].label.second = root;
    }

    ::tree<Problem> build_tree(std::vector<int> component, std::map<int, std::vector<Edge> >& emb) {
        int level = vertex_level[component[0]];

//        std::reverse(component.begin(), component.end());

        std::map<graph_traits<Graph>::edge_descriptor, std::vector<int> > faces;
        std::vector<std::vector<int> > vertices_in_face;
        face_getter<Edge> my_vis(&faces, vertices_in_face);
        level_face_traversal<Graph>(emb, my_vis);
        ::tree<Problem> t(my_vis.current_face, level);

        int outer = 0;
        for (; outer < vertices_in_face.size(); outer++) {
            auto& face = vertices_in_face[outer];

            int start = -1;
            for (int i = 0; i < face.size(); i++) {
                if (face[i] == component[0]) {
                    start = i;
                    break;
                }
            }
            if (start == -1) {
                continue;
            }

            std::rotate(face.begin(), face.begin() + start, face.end());

            if (face == component) {
                break;
            }
        }

        t.outer_face = outer;
        tree_builder<graph_traits<Graph>::edge_descriptor, Problem, PlanarEmbedding> tree_b(faces, t, g, embedding);
        level_face_traversal<Graph>(emb, tree_b);
        t.remove_outer_face();

        for (auto& node : t.t) {
            std::reverse(node.children.begin(), node.children.end());
            std::reverse(node.face.begin(), node.face.end());
        }

        return t;
    }

    ::tree<Problem> build_tree_with_dividing_points(std::vector<int> component, int root) {
        int level = vertex_level[component[0]];

        if (component.size() == 1) {
            ::tree<Problem> t(1, level);
            t[0].label.first = component[0];
            t[0].label.second = component[0];
            t.root = 0;
            return t;
        }

        Graph& g_temp = g.create_subgraph(component.begin(), component.end());

        auto bicomp = get(edge_index, g_temp);
        int bi_num = biconnected_components(g_temp, bicomp);
        std::vector< ::tree<Problem> > trees(bi_num);
        std::vector< std::set<int> > v_in_bicomps(bi_num);

        std::map<std::pair<int, int>, int> bi_map;
        graph_traits<Graph>::edge_iterator ei, ei_end;
        for(boost::tie(ei, ei_end) = edges(g_temp); ei != ei_end; ++ei) {
            int global_source = g_temp.local_to_global(ei->m_source);
            int global_target = g_temp.local_to_global(ei->m_target);
            v_in_bicomps[bicomp[*ei]].insert(global_source);
            v_in_bicomps[bicomp[*ei]].insert(global_target);
            bi_map[std::pair<int, int>(global_source, global_target)] = bicomp[*ei];
            bi_map[std::pair<int, int>(global_target, global_source)] = bicomp[*ei];
        }

        std::vector< std::vector<int> > bicomps(bi_num);
        for (int i = 0; i < component.size(); i++) {
            int curr = component[i];
            int next = component[(i + 1) % component.size()];
            bicomps[bi_map[std::pair<int, int>(curr, next)]].push_back(curr);
        }

        for (auto& bic : bicomps) {
            if (bic.size() == 2) {
                Edge e = add_edge(bic[0], bic[1], g).first;
                int it = get_edge_it(e, bic[0], embedding);
                embedding[bic[0]].insert(embedding[bic[0]].begin() + it, e);
                it = get_edge_it(e, bic[1], embedding);
                embedding[bic[1]].insert(embedding[bic[1]].begin() + it, e);
            }
        }

        for (int i = 0; i < bi_num; i++) {
            std::map<int, std::vector<Edge> > emb;
            for (int v : v_in_bicomps[i]) {
                for (Edge e : embedding[v]) {
                    if (bi_map[std::pair<int, int>(e.m_source, e.m_target)] == i &&
                    vertex_level[e.m_source] == level && vertex_level[e.m_target] == level) {
                        emb[v].push_back(e);
                    }
                }
            }

            trees[i] = build_tree(bicomps[i], emb);
        }

        std::vector<int> art_points;
        articulation_points(g_temp, std::back_inserter(art_points));

        for (int& a : art_points) {
            a = g_temp.local_to_global(a);
        }

        std::map<int, std::vector<int> > art_to_bicomp;
        std::map<int, std::vector<int> > bicomp_to_art;

        for (int a : art_points) {
            for (int i = 0; i < v_in_bicomps.size(); i++) {
                if (v_in_bicomps[i].find(a) != v_in_bicomps[i].end()) {
                    art_to_bicomp[a].push_back(i);
                    bicomp_to_art[i].push_back(a);
                }
            }
        }

        int main_bicomp = 0;
        for (int i = 0; i < v_in_bicomps.size(); i++) {
            if (v_in_bicomps[i].find(root) != v_in_bicomps[i].end()) {
                main_bicomp = i;
                break;
            }
        }

        root_tree_with_root(trees[main_bicomp], root);

        std::vector<bool> merged(bi_num, false);
        merged[main_bicomp] = true;

        std::vector<int> where_to_merge(bi_num);
        where_to_merge[main_bicomp] = trees[main_bicomp].root;

        std::queue<int> q;
        q.push(main_bicomp);

        std::map<int, int> place_in_comp;

        for (int i = 0; i < component.size(); i++) {
            if (place_in_comp.find(component[i]) == place_in_comp.end()) {
                place_in_comp[component[i]] = i;
            }
        }

        while(!q.empty()) {
            int cur_bicomp = q.front();
            q.pop();

            for (int a : bicomp_to_art[cur_bicomp]) {
                for (int b : art_to_bicomp[a]) {
                    if (!merged[b]) {
                        q.push(b);
                        root_tree_with_root(trees[b], a);

                        merged[b] = true;
                        where_to_merge[b] = trees[main_bicomp].merge(trees[b], place_in_comp,
                                                                     where_to_merge[cur_bicomp]);
                    }
                }
            }
        }


        for (int i = 0; i < trees[main_bicomp].size(); i++) {
            auto& t = trees[main_bicomp][i];

            if (t.children.empty()) {
                continue;
            }

            auto& face = t.face;

            std::vector< std::pair<int, int> > v_in_c;
            check_for_components(face, v_in_c);
            if (!v_in_c.empty()) {
                std::vector< std::vector<int> > components;
                std::map<int, int> v_to_c;
                get_component(components, v_to_c, v_in_c);

                if (components.size() > 1) {
                    make_connected(components, v_to_c, v_in_c);
                    components.clear();
                    v_to_c.clear();
                    get_component(components, v_to_c, v_in_c);
                }

                int face_v = v_in_c[0].second;
                for (auto e : v_in_c) {
                    if (e.second != face_v) {
                        face_v = -1;
                        break;
                    }
                }

                if (face_v > -1) {
                    int next_in_face;
                    for (int v = 0; v < face.size(); v++) {
                        if (face[v] == face_v) {
                            face_v = face[(v + 1) % face.size()];
                            next_in_face = face[(v + 2) % face.size()];
                            break;
                        }
                    }

                    int comp_v = v_in_c.back().first;
                    Edge new_e = add_edge(face_v, comp_v, g).first;
                    added_edges.emplace(face_v, comp_v);
                    embedding[face_v].insert(embedding[face_v].begin() + get_edge_it(face_v, next_in_face, embedding) + 1, new_e);
                    embedding[comp_v].insert(embedding[comp_v].begin() + get_edge_it(comp_v, face_v, embedding) + 1, new_e);
                }

                triangulate(face, components[0], 1);

                std::rotate(face.begin(), std::find(face.begin(), face.end(), t.label.first), face.end());
                v_in_c.clear();
                check_for_components(face, v_in_c);
                components.clear();
                v_to_c.clear();
                components.emplace_back();
//                get_component(components, v_to_c, v_in_c);


                int v;
                if (t.label.first != t.label.second) {
                    v = find_third(t.label.first, t.label.second, components[0]);
                } else {
                    int w;
                    for (int c = t.children.size() - 1; c >= 0; c--) {
                        if (trees[main_bicomp][t.children[c]].label.first != t.label.first) {
                            w = trees[main_bicomp][t.children[c]].label.first;
                            break;
                        }
                    }
                    v = find_third(t.label.first, w, components[0]);
                }
                t.component_tree = build_tree_with_dividing_points(components[0], v);
                t.component_tree.enclosing_tree = &trees[main_bicomp];
                t.component_tree.enclosing_face = i;
            }
        }

        return trees[main_bicomp];
    }

    void create_boudaries_rec(::tree<Problem>& t, int node) {
        if (t[node].children.empty()) {
            return;
        }

        for (int i : t[node].children) {
            create_boudaries_rec(t, i);
        }

        t[node].LB = t[t[node].children[0]].LB;
        t[node].RB = t[t[node].children.back()].RB;
    }

    void create_boudaries(::tree<Problem>& t, ::tree<Problem>& t2, int root) {
        if (t.enclosing_tree != nullptr) {
            Problem& f = t.get_enclosing_face();
            std::vector<int> y_table;
            y_table.push_back(t2[f.children[0]].label.first);

            for (int v : f.children) {
                y_table.push_back(t2[v].label.second);
            }

            std::vector<int> leaves;
            get_leaves(t, leaves, root);
            t[leaves[0]].LB = 0;
            t[leaves.back()].RB = f.children.size();

            for (int j = 1; j < leaves.size(); j++) {
                Problem &v = t[leaves[j]];
                Problem &w = t[leaves[j - 1]];
                Edge one(w.label.first, w.label.second, nullptr);
                Edge two(v.label.first, v.label.second, nullptr);

                std::set<int> dividing_points;
                int div = find_dividing_points(one, two, dividing_points);

                for (int i = w.LB; i < y_table.size(); i++) {
                    if (y_table[i] == div) {
                        v.LB = w.RB = i;
                        break;
                    }
                }
            }

            create_boudaries_rec(t, root);
        } else {
            for (auto &node : t.t) {
                node.LB = node.label.first;
                node.RB = node.label.second;
            }
        }

        for (int i = 0; i < t.size(); i++) {
            Problem& node = t[i];
            if (!node.component_tree.empty()) {
                create_boudaries(node.component_tree, t, node.component_tree.root);
            }
        }
    }

    void table (::tree<Problem>& t, int v) {
        int level = t[v].level;

        if (!t[v].children.empty() && t[v].component_tree.empty()) {
            table(t, t[v].children[0]);
            t[v].val = t[t[v].children[0]].val;

            for(int i=1; i < t[v].children.size(); i++) {
                table(t, t[v].children[i]);
                t[v].merge(t[v].val, t[t[v].children[i]].val);
            }

            t[v].adjust(g, added_edges);
            return;
        } else if (!t[v].children.empty() && !t[v].component_tree.empty()) {
            auto& ct = t[v].component_tree;
            table(ct, ct.root);
            t[v].contract(ct[ct.root], g, added_edges);
            t[v].adjust(g, added_edges);
            return;
        } else if (level > 1) {
            std::vector<int> z_table;
            z_table.push_back(t.get_enclosing_face().get_child(0).label.first);
            for (int x : t.enclosing_tree->t[t.enclosing_face].children) {
                z_table.push_back(t.enclosing_tree->t[x].label.second);
            }

            if (z_table.front() == z_table.back()) {
                z_table.pop_back();
            }

            Edge& e = embedding[t[v].label.first][get_edge_it(t[v].label.first, t[v].label.second, embedding) - 1];
            int third = e.m_source == t[v].label.first ?  e.m_target : e.m_source;
            int p = t[v].LB;
            for (int i = t[v].RB; i > t[v].LB; i--) {
                if (i >= z_table.size()) {
                    i--;
                }
                if (z_table[i] == third) {
                    p = i;
                    break;
                }
            }

            t[v].create(p, g, added_edges);
            t[v].adjust(g, added_edges);

            int j = p - 1;

            while (j >= t[v].LB) {
                table(*t.enclosing_tree, t.enclosing_tree->t[t.enclosing_face].children[j]);
                Problem second = t.get_enclosing_face().get_child(j).extend(t[v].label.first, g, added_edges);
                t[v].merge(second.val, t[v].val);
                j--;
            }

            j = p;

            while (j < t[v].RB) {
                table(*t.enclosing_tree, t.enclosing_tree->t[t.enclosing_face].children[j]);
                Problem second = t.get_enclosing_face().get_child(j).extend(t[v].label.second, g, added_edges);
                t[v].merge(t[v].val, second.val);
                j++;
            }
        }
    }

public:
    baker_impl(const Graph& arg_g, PlanarEmbedding emb, std::vector<int> out_face): g(arg_g), embedding(emb),
    outer_face(out_face), vertex_level(num_vertices(arg_g)) {
        std::vector< std::vector<Edge> > outer_edges;
        int k = name_levels(embedding, outer_face,vertex_level, outer_edges);

        int v = 0;

        for (int i = 0; i < vertex_level.size(); i++) {
            if (vertex_level[i] == 1) {
                v = i;
                break;
            }
        }

        tree = build_tree_with_dividing_points(outer_face, v);

        create_boudaries(tree, tree, tree.root);

        table(tree, tree.root);
    }

    int result() {
        return tree[tree.root].result();
    }

};

template <typename Problem>
int baker(const Graph& g, PlanarEmbedding embedding, std::vector<int> outer_face) {
    baker_impl < Problem > b(g, embedding, outer_face);
    return b.result();
}

#define TECHNIKA_BAKER_BAKER_K_OUTER_PLANAR_HPP

#endif //TECHNIKA_BAKER_BAKER_K_OUTER_PLANAR_HPP
