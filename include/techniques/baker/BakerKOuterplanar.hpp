/*
 * BakerKOuterplanar.hpp
 *
 *  Created on: 08.02.2021
 *      Author: Mikołaj Twaróg
 *      Ported by: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

#include <iostream>
#include <vector>
#include <climits>
#include <queue>

#include <boost/graph/biconnected_components.hpp>

#include <structures/CyclicVector.hpp>
#include <techniques/baker/LevelFaceTraversal.hpp>
#include <techniques/baker/NameLevels.hpp>
#include <techniques/baker/PlanarTree.hpp>
#include <techniques/baker/Visitors.hpp>

template<typename Problem>
class baker_impl {
 public:
    baker_impl(const Graph& arg_g, const NetworKit::Graph &graph, PlanarEmbedding &emb, const std::vector<int> &out_face)
            : graph(std::make_optional(graph)), g(arg_g), embedding(convert(emb)), outer_face(out_face), vertex_level(boost::num_vertices(arg_g)) {
        auto [k, _] = name_levels(embedding, outer_face, vertex_level);
        int v = std::find_if(
            vertex_level.begin(), vertex_level.end(),
            [](int level) { return level == 1; }) - vertex_level.begin();
        tree = build_tree_with_dividing_points(outer_face, v);
        create_boundaries(tree, tree, tree.root);
        table(tree, tree.root);
    }

    int result() {
        return tree[tree.root].result();
    }

 private:
    std::optional<NetworKit::Graph> graph;

    Graph g;
    PlanarEmbedding2 embedding;
    std::vector<int> outer_face, vertex_level;
    PlanarTree<Problem> tree;
    std::set<std::pair<int, int>> added_edges;

    void triangulate(std::vector<int> &face, std::vector<int> &component) {
        if (component.size() == 1) {
            int c = component[0];
            embedding[c].clear();
            for (int i = 0; i < face.size(); i++) {
                int v = face[i];
                embedding[c].emplace_back(v, c);
                if (!edge(v, c, g).second) {
                    Edge new_e = add_edge(c, v, g).first;
                    auto new_e2 = NetworKit::Edge(new_e.m_source, new_e.m_target);
                    int next = face[(i + 1) % face.size()];
                    embedding[v].insert(
                        embedding[v].begin() + get_edge_it(v, next, embedding) + 1, new_e2);
                    added_edges.emplace(v, c);
                }
            }
            return;
        }

        int level = vertex_level[face[0]];
        std::optional<NetworKit::Edge> connecting_e;
        for (int i = 0; i < face.size(); i++) {
            int v = face[i];
            for (auto &e : embedding[v]) {
                int neighbour = e.u == v ? e.v : e.u;
                auto element = std::find(component.begin(), component.end(), neighbour);
                if (element != component.end()) {
                    connecting_e = e;
                    std::rotate(face.begin(), face.begin() + i, face.end());
                    std::rotate(component.begin(), element, component.end());
                    break;
                }
            }
            if (connecting_e) {
                break;
            }
        }
        int face_it = 0, comp_it = 0;
        int curr_e_it_face = get_edge_it(*connecting_e, face.front(), embedding);
        int curr_e_it_comp = get_edge_it(*connecting_e, component.front(), embedding);
        std::set<int> comp_vis, face_vis;
        int comp = 0;
        while (comp < component.size() || face_vis.size() < face.size()) {
            bool res = false;
            int comp_curr = component[comp_it], face_curr = face[face_it], last = -1;

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
                auto &e = embedding[comp_curr][i];
                int neighbour = e.u == comp_curr ? e.v : e.u;
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
                    NetworKit::Edge e(face[i - 1], face[i]);

                    int e_it;
                    if (i == face_it) {
                        e_it = curr_e_it_face;
                    } else {
                        e_it = get_edge_it(e, face[i], embedding);
                    }


                    if (embedding[face[i]][e_it].u != comp_curr
                            && embedding[face[i]][e_it].v != comp_curr) {
                        auto new_e = add_edge(face[i], comp_curr, g).first;
                        auto new_e2 = NetworKit::Edge(new_e.m_source, new_e.m_target);
                        embedding[face[i]].insert(embedding[face[i]].begin() + e_it, new_e2);
                        embedding[comp_curr].insert(embedding[comp_curr].begin() + curr_e_it_comp, new_e2);
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

                if (embedding[comp_curr][temp_e_it_comp].u != face_curr && embedding[comp_curr][temp_e_it_comp].v != face_curr) {
                    auto new_e = add_edge(face_curr, comp_curr, g).first;
                    auto new_e2 = NetworKit::Edge(new_e.m_source, new_e.m_target);
                    embedding[face_curr].insert(
                        embedding[face_curr].begin() + curr_e_it_face, new_e2);
                    embedding[comp_curr].insert(
                        embedding[comp_curr].begin() + curr_e_it_comp + 1, new_e2);
                } else {
                    auto &e_con = embedding[comp_curr][temp_e_it_comp];
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

    int find_third(int one, int two, std::vector<int> &component) {
        auto &one_edges = embedding[one];
        int level = vertex_level[one];
        if (one == two) {
            for (auto e_i : one_edges) {
                int target_i = e_i.u == one ? e_i.v : e_i.u;
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
            third = one_edges[i].u == one ? one_edges[i].v : one_edges[i].u;

            if (third != one) {
                for (int j = 0; j < embedding[third].size(); j++) {
                    auto curr_e = embedding[third][j];
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
        NetworKit::Edge current_e;
        int current_v;
        bool res = false;
        level++;
        auto &edges = embedding[starting_v];
        for (int i = (connecting_e_it + 1) % edges.size(); i != connecting_e_it; i = (i + 1) % edges.size()) {
            auto e = edges[i];
            int neighbour = e.u == starting_v ? e.v : e.u;
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
            res = false;
            for (int i = (current_e_it + 1 + edges2.size()) % edges2.size(); i != current_e_it; i = (i + 1 + edges2.size()) % edges2.size()) {
                auto e = edges2[i];
                int neighbour = e.u == current_v ? e.v : e.u;

                if (vertex_level[neighbour] == level) {
                    current_e = e;
                    current_v = neighbour;
                    res = true;
                    break;
                }
            }
            if (res == false) {
                auto e = edges2[current_e_it];
                int neighbour = e.u == current_v ? e.v : e.u;

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

    int find_dividing_points(NetworKit::Edge one, NetworKit::Edge two) {
        std::set<int> dividing_points;
        auto v = one.v;
        auto level = vertex_level[v];
        auto one_it = get_edge_it(one, v, embedding), two_it = get_edge_it(two, v, embedding);
        auto &edges = embedding[v];
        for (int i = (one_it + 1) % edges.size(); i != two_it; i = (i + 1) % edges.size()) {
            auto e = edges[i];
            if (vertex_level[e.u] == level - 1 || vertex_level[e.v] == level - 1) {
                 dividing_points.insert(vertex_level[e.u] == level - 1 ? e.u : e.v);
                 return vertex_level[e.u] == level - 1 ? e.u : e.v;
            }
        }
        return -1;
    }

    void get_leaves(PlanarTree<Problem>& t, std::vector<int> &leaves, int v) {
        if (t[v].children.empty()) {
            leaves.push_back(v);
            return;
        }
        for (auto child : t[v].children) {
            get_leaves(t, leaves, child);
        }
    }

    void get_component(
            std::vector<std::vector<int>> &components, std::map<int, int> &vis,
            cyclic_vector<std::pair<int, int>> &v_in_c) {
        if (vertex_level[v_in_c[0].first] == 1) {
            components[0].insert(components[0].end(), outer_face.begin(), outer_face.end());
            return;
        }
        for (int c = 0, comp_num = 0; c < v_in_c.size(); c++) {
            if (vis.find(v_in_c[c].first) != vis.end()) {
                continue;
            }
            int starting_v = v_in_c[c].first;
            int level = vertex_level[starting_v];
            int connecting_e_it = get_edge_it(starting_v, v_in_c[c].second, embedding);

            vis[starting_v] = comp_num;
            components.emplace_back();
            components.back().push_back(starting_v);

            NetworKit::Edge current_e;
            int current_v;
            bool res = false;
            auto &edges = embedding[starting_v];
            for (int i = (connecting_e_it + 1) % edges.size(); i != connecting_e_it; i = (i + 1) % edges.size()) {
                auto e = edges[i];
                int neighbour = e.u == starting_v ? e.v : e.u;
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
                res = false;
                for (int i = (current_e_it + 1 + edges2.size()) % edges2.size(); i != current_e_it;
                        i = (i + 1 + edges2.size()) % edges2.size()) {
                    auto e = edges2[i];
                    int neighbour = e.u == current_v ? e.v : e.u;

                    if (vertex_level[neighbour] == level) {
                        current_e = e;
                        current_v = neighbour;
                        res = true;
                        break;
                    }
                }
                if (res == false) {
                    auto e = edges2[current_e_it];
                    int neighbour = e.u == current_v ? e.v : e.u;
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

    void make_connected(
            std::vector<std::vector<int>> &components, std::map<int, int> &vis,
            cyclic_vector<std::pair<int, int>> &v_in_c) {
        std::vector<bool> connected(components.size(), false);
        for (int i = 0; i < v_in_c.size(); i++) {
            auto [curr, curr_con] = v_in_c[i];
            auto [next, next_con] = v_in_c[i + 1];
            connected[vis[curr]] = true;
            if (vis[curr] == vis[next] || edge(curr, next, g).second || connected[vis[next]]) {
                continue;
            }
            added_edges.emplace(curr, next);
            Edge new_e = add_edge(curr, next, g).first;
            auto new_e2 = NetworKit::Edge(new_e.m_source, new_e.m_target);
            embedding[curr].insert(embedding[curr].begin() + get_edge_it(curr, curr_con, embedding) + 1, new_e2);
            embedding[next].insert(embedding[next].begin() + get_edge_it(next, next_con, embedding), new_e2);
        }
    }

    void check_for_components(const std::vector<int> &face, std::vector<std::pair<int, int>>& out) {
        if (face.size() == 2) {
            return;
        }
        int level = vertex_level[face[0]], prev = face.size() - 1, curr = 0, next = 1;
        for (; curr < face.size(); prev = (prev + 1) % face.size(), curr++, next = (next + 1) % face.size()) {
            int one = get_edge_it(face[curr], face[prev], embedding);
            int two = get_edge_it(face[curr], face[next], embedding);
            for (int i = (one - 1 + embedding[face[curr]].size()) % embedding[face[curr]].size(); i != two;
                    i = (i - 1 + embedding[face[curr]].size()) % embedding[face[curr]].size()) {
                auto e = embedding[face[curr]][i];
                int target = e.u == face[curr] ? e.v : e.u;
                if (vertex_level[target] == level + 1) {
                    out.emplace_back(target, face[curr]);
                }
            }
        }
    }

    void root_tree_2(PlanarTree<Problem> &t, int node) {
        if (t[node].children.empty()) {
            return;
        }
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
        std::rotate(
            t[node].children.begin(), t[node].children.begin() + parent_it, t[node].children.end());
        t[node].label.first = t[t[node].children[0]].label.first;
        t[node].label.second = t[t[node].children.back()].label.second;
    }

    void root_tree_with_root(PlanarTree<Problem> &t, int root) {
        int node = std::find_if(
            t.t.begin(), t.t.end(),
            [&] (auto v) { return std::find(v.face.begin(), v.face.end(), root) != v.face.end();
        }) - t.t.begin();
        t.root = node;
        root_tree_2(t, node);
        auto element = std::find_if(
            t[node].children.begin(), t[node].children.end(),
            [&] (auto child) { return t[child].label.first == root; });
        std::rotate(t[node].children.begin(), element, t[node].children.end());
        t[node].label.first = t[node].label.second = root;
    }

    PlanarTree<Problem> build_tree(std::vector<int> component, PlanarEmbedding2 &embedding) {
        int level = vertex_level[component[0]];
        std::map<std::pair<int, int>, std::vector<int>> faces;
        std::vector<std::vector<int>> vertices_in_face;
        face_getter visitor(faces, vertices_in_face);
        level_face_traversal(embedding, visitor);
        PlanarTree<Problem> t(visitor.current_face, level);

        int outer = 0;
        for (; outer < vertices_in_face.size(); outer++) {
            auto &face = vertices_in_face[outer];
            auto element = std::find_if(
                face.begin(), face.end(), [&] (auto c) { return c == component[0]; });
            if (element == face.end()) {
                continue;
            }
            std::rotate(face.begin(), element, face.end());
            if (face == component) {
                break;
            }
        }
        t.outer_face = outer;
        auto tree_b = tree_builder<Problem>(faces, t, g);
        level_face_traversal(embedding, tree_b);
        t.remove_outer_face();

        for (auto &node : t.t) {
            std::reverse(node.children.begin(), node.children.end());
            std::reverse(node.face.begin(), node.face.end());
        }
        return t;
    }

    PlanarTree<Problem> build_tree_with_dividing_points(
            const std::vector<int> &component, int root) {
        int level = vertex_level[component[0]];

        if (component.size() == 1) {
            PlanarTree<Problem> t(1, level);
            t[0].label.first = t[0].label.second = component[0];
            t.root = 0;
            return t;
        }

        Graph& g_temp = g.create_subgraph(component.begin(), component.end());
        auto bicomp = get(boost::edge_index, g_temp);
        int bi_num = biconnected_components(g_temp, bicomp);
        std::vector<PlanarTree<Problem>> trees(bi_num);
        std::vector<std::set<int>> v_in_bicomps(bi_num);
        std::unordered_map<NetworKit::Edge, int> bi_map;
        boost::graph_traits<Graph>::edge_iterator ei, ei_end;
        for(boost::tie(ei, ei_end) = edges(g_temp); ei != ei_end; ++ei) {
            int global_source = g_temp.local_to_global(ei->m_source);
            int global_target = g_temp.local_to_global(ei->m_target);
            v_in_bicomps[bicomp[*ei]].insert(global_source);
            v_in_bicomps[bicomp[*ei]].insert(global_target);
            bi_map[NetworKit::Edge(global_source, global_target)] = bicomp[*ei];
            bi_map[NetworKit::Edge(global_target, global_source)] = bicomp[*ei];
        }

        std::vector<std::vector<int>> bicomps(bi_num);
        for (int i = 0; i < component.size(); i++) {
            int curr = component[i], next = component[(i + 1) % component.size()];
            bicomps[bi_map[NetworKit::Edge(curr, next)]].push_back(curr);
        }
        for (auto &bic : bicomps) {
            if (bic.size() == 2) {
                auto e = add_edge(bic[0], bic[1], g).first;
                auto e2 = NetworKit::Edge(e.m_source, e.m_target);
                for (auto v : {bic[0], bic[1]}) {
                    auto pos = get_edge_it(e2, v, embedding);
                    if (embedding[v][pos].u == e2.u) {
                        std::swap(e2.u, e2.v);
                    }
                    embedding[v].insert(embedding[v].begin() + pos, e2);
                }
            }
        }

        for (int i = 0; i < bi_num; i++) {
            PlanarEmbedding2 emb(graph->numberOfNodes());
            for (int v : v_in_bicomps[i]) {
                for (const auto &e : embedding[v]) {
                    if (bi_map[e] == i && vertex_level[e.u] == level && vertex_level[e.v] == level) {
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

        std::map<int, std::vector<int>> art_to_bicomp, bicomp_to_art;
        for (int a : art_points) {
            for (int i = 0; i < v_in_bicomps.size(); i++) {
                if (v_in_bicomps[i].find(a) != v_in_bicomps[i].end()) {
                    art_to_bicomp[a].push_back(i);
                    bicomp_to_art[i].push_back(a);
                }
            }
        }
        int main_bicomp = std::find_if(
            v_in_bicomps.begin(), v_in_bicomps.end(),
            [&] (auto v) { return v.find(root) != v.end(); }) - v_in_bicomps.begin();

        root_tree_with_root(trees[main_bicomp], root);

        std::vector<bool> merged(bi_num, false);
        merged[main_bicomp] = true;

        std::vector<int> where_to_merge(bi_num);
        where_to_merge[main_bicomp] = trees[main_bicomp].root;

        std::queue<int> q({main_bicomp});
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
                        where_to_merge[b] = trees[main_bicomp].merge(
                            trees[b], place_in_comp, where_to_merge[cur_bicomp]);
                    }
                }
            }
        }
        for (int i = 0; i < trees[main_bicomp].size(); i++) {
            auto &t = trees[main_bicomp][i];
            if (t.children.empty()) {
                continue;
            }
            auto &face = t.face;
            cyclic_vector<std::pair<int, int>> v_in_c;
            check_for_components(face, v_in_c);
            if (!v_in_c.empty()) {
                std::vector<std::vector<int>> components;
                std::map<int, int> v_to_c;
                get_component(components, v_to_c, v_in_c);
                if (components.size() > 1) {
                    make_connected(components, v_to_c, v_in_c);
                    components.clear(), v_to_c.clear();
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
                    auto new_e2 = NetworKit::Edge(new_e.m_source, new_e.m_target);
                    embedding[face_v].insert(embedding[face_v].begin() + get_edge_it(face_v, next_in_face, embedding) + 1, new_e2);
                    embedding[comp_v].insert(embedding[comp_v].begin() + get_edge_it(comp_v, face_v, embedding) + 1, new_e2);
                }

                triangulate(face, components[0]);

                std::rotate(face.begin(), std::find(face.begin(), face.end(), t.label.first), face.end());
                v_in_c.clear();
                check_for_components(face, v_in_c);
                components.clear();
                components.emplace_back();
                int v;
                if (t.label.first != t.label.second) {
                    v = find_third(t.label.first, t.label.second, components[0]);
                } else {
                    int w2;
                    for (int c = t.children.size() - 1; c >= 0; c--) {
                        if (trees[main_bicomp][t.children[c]].label.first != t.label.first) {
                            w2 = trees[main_bicomp][t.children[c]].label.first;
                            break;
                        }
                    }
                    auto element = std::find_if(
                        t.children.rbegin(), t.children.rend(),
                        [&](auto c) {
                            return trees[main_bicomp][c].label.first != t.label.first;
                        });
                    assert(w2 == trees[main_bicomp][*element].label.first);
                    v = find_third(t.label.first, w2, components[0]);
                }
                t.component_tree = build_tree_with_dividing_points(components[0], v);
                t.component_tree.enclosing_tree = &trees[main_bicomp];
                t.component_tree.enclosing_face = i;
            }
        }
        return trees[main_bicomp];
    }

    void create_boundaries_rec(PlanarTree<Problem>& t, int node) {
        if (t[node].children.empty()) {
            return;
        }
        for (auto child : t[node].children) {
            create_boundaries_rec(t, child);
        }
        t[node].left = t[t[node].children[0]].left;
        t[node].right = t[t[node].children.back()].right;
    }

    void create_boundaries(PlanarTree<Problem> &t, PlanarTree<Problem> &t2, int root) {
        if (t.enclosing_tree != nullptr) {
            auto &f = t.get_enclosing_face();
            std::vector<int> y_table{t2[f.children[0]].label.first};
            for (int v : f.children) {
                y_table.push_back(t2[v].label.second);
            }
            std::vector<int> leaves;
            get_leaves(t, leaves, root);
            t[leaves[0]].left = 0, t[leaves.back()].right = f.children.size();
            for (int j = 1; j < leaves.size(); j++) {
                auto &v = t[leaves[j]], &w = t[leaves[j - 1]];
                NetworKit::Edge one(w.label.first, w.label.second);
                NetworKit::Edge two(v.label.first, v.label.second);
                int div = find_dividing_points(one, two);
                for (int i = w.left; i < y_table.size(); i++) {
                    if (y_table[i] == div) {
                        v.left = w.right = i;
                        break;
                    }
                }
            }
            create_boundaries_rec(t, root);
        } else {
            for (auto &node : t.t) {
                node.left = node.label.first, node.right = node.label.second;
            }
        }
        for (auto &node : t.t) {
            if (!node.component_tree.empty()) {
                create_boundaries(node.component_tree, t, node.component_tree.root);
            }
        }
    }

    void table(PlanarTree<Problem>& t, int v) {
        int level = t[v].level;
        if (!t[v].children.empty() && t[v].component_tree.empty()) {
            table(t, t[v].children[0]);
            t[v].val = t[t[v].children[0]].val;
            for(int i = 1; i < t[v].children.size(); i++) {
                table(t, t[v].children[i]);
                t[v].merge(t[v].val, t[t[v].children[i]].val);
            }
            t[v].adjust(g, added_edges);
            return;
        }
        if (!t[v].children.empty() && !t[v].component_tree.empty()) {
            auto &ct = t[v].component_tree;
            table(ct, ct.root);
            t[v].contract(ct[ct.root], g, added_edges);
            t[v].adjust(g, added_edges);
            return;
        }
        if (level > 1) {
            std::vector<int> z_table;
            z_table.push_back(t.get_enclosing_face().get_child(0).label.first);
            for (int x : t.get_enclosing_face().children) {
                z_table.push_back(t.enclosing_tree->t[x].label.second);
            }
            if (z_table.front() == z_table.back()) {
                z_table.pop_back();
            }

            auto &e = embedding[t[v].label.first][get_edge_it(t[v].label.first, t[v].label.second, embedding) - 1];
            int third = e.u == t[v].label.first ? e.v : e.u;
            int p = t[v].left;
            for (int i = t[v].right; i > t[v].left; i--) {
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
            for (int j = p - 1; j >= t[v].left; j--) {
                table(*t.enclosing_tree, t.get_enclosing_face().children[j]);
                Problem second = t.get_enclosing_face().get_child(j).extend(t[v].label.first, g, added_edges);
                t[v].merge(second.val, t[v].val);
            }
            for (int j = p; j < t[v].right; j++) {
                table(*t.enclosing_tree, t.get_enclosing_face().children[j]);
                Problem second = t.get_enclosing_face().get_child(j).extend(t[v].label.second, g, added_edges);
                t[v].merge(t[v].val, second.val);
            }
        }
    }
};

template <typename Problem>
int baker(const Graph& g, const NetworKit::Graph &graph, PlanarEmbedding embedding, std::vector<int> outer_face) {
    baker_impl<Problem> b(g, graph, embedding, outer_face);
    return b.result();
}
