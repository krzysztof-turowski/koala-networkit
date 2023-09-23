/*
 * CreateTreeDecomposition.hpp
 *
 *  Created on: 18.05.2021
 *      Author: Mikołaj Twaróg
 *      Ported by: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

#include <list>
#include <map>
#include <queue>
#include <set>
#include <vector>

#include <techniques/baker/LevelFaceTraversal.hpp>
#include <techniques/baker/NameLevels.hpp>
#include <techniques/baker/Visitors.hpp>

struct tree_decomposition {
    std::vector<cyclic_vector<int>> tree;
    std::vector<std::set<int>> nodes;
    int n;

    std::map<std::pair<int, int>, int> edge_it;

    std::set<int>& operator[](int x) {
        return nodes[x];
    }

    std::set<int>& operator[](std::pair<int, int> e) {
        return nodes[edge_it[e]];
    }

    std::set<int>& operator[](NetworKit::Edge e) {
        return nodes[edge_it[std::minmax(e.u, e.v)]];
    }

    std::set<int>& operator[](Edge e) {
        return nodes[edge_it[std::minmax(e.m_source, e.m_target)]];
    }
};

class create_tree_decomposition {
    Graph graph;
    std::vector<cyclic_vector<Edge>> embedding;
    std::vector<cyclic_vector<NetworKit::Edge>> embedding2;
    std::vector<int> outer_face;
    tree_decomposition tr;
    std::vector<int> vertex_level;
    std::vector<cyclic_vector<int>> spanning_forest;
    int k;
    std::vector<std::vector<cyclic_vector<int>>> level_graphs, level_cycles;
    int outer_face_it;
    std::vector<int> expanded_vertices, grey;
    int cycle_num = 0;

    int cycle_dfs(int v, int p, int level, std::vector<int> &visited) {
        visited[v] = 1;
        int cycle = -1;
        for (int nei : level_graphs[level][v]) {
            if (visited[nei] == 1 && nei != p) {
                grey.push_back(nei);
                cycle = cycle_num++;
                level_cycles[level].emplace_back();
            }
            if (visited[nei] == 0) {
                cycle = std::max(cycle, cycle_dfs(nei, v, level, visited));
            }
        }
        visited[v] = 2;
        if (cycle >= 0) {
            level_cycles[level][cycle].push_back(v);
        }
        if (!grey.empty() && grey.back() == v) {
            cycle = -1;
            grey.pop_back();
        }
        return cycle;
    }

    bool check_for_edge(int v, int w, std::vector<cyclic_vector<int>>& g) {
        return std::any_of(g[v].begin(), g[v].end(), [w](auto x) { return x == w; });
    }

    bool check_for_edge(int v, int w, std::vector<cyclic_vector<Edge>>& g) {
        for (Edge& e : g[v]) {
            int nei = v == e.m_source ? e.m_target : e.m_source;
            if (nei == w) {
                return true;
            }
        }
        return false;
    }

    void get_leaves(int v, int p, int level, std::set<int>& leaves) {
        for (int w : spanning_forest[v]) {
            if (w != p) {
                get_leaves(w, v, level, leaves);
            }
        }
        if (vertex_level[v] == level) {
            leaves.insert(v);
        }
    }

    void outerplanar() {
        std::vector<std::list<int>> g(embedding.size());
        std::stack<int> low_degree;
        for (int i = 0; i < embedding.size(); i++) {
            if (embedding[i].size() <= 2) {
                low_degree.push(i);
            }

            for (Edge& edge : embedding[i]) {
                int v = edge.m_source == i ? edge.m_target : edge.m_source;
                g[i].push_back(v);
            }
        }

        std::vector<std::vector<int>> deleted;
        for (int i = 0; i < g.size() - 2; i++) {
            int v = low_degree.top();
            low_degree.pop();

            if (g[v].size() == 0) {
                i--;
                continue;
            }

            deleted.emplace_back();
            deleted.back().push_back(v);

            if (g[v].size() == 1) {
                int nei = g[v].front();
                deleted.back().push_back(nei);

                for (auto it = g[nei].begin(); it != g[nei].end(); ) {
                    if (*it == v) {
                        it = g[nei].erase(it);
                    } else {
                        ++it;
                    }
                }
                if (g[nei].size() == 2) {
                    low_degree.push(nei);
                }
                g[v].clear();
                continue;
            }

            int n_one = g[v].front(), n_two = g[v].back();
            deleted.back().push_back(n_one);
            deleted.back().push_back(n_two);
            bool edge_needed = true;

            for (auto it = g[n_one].begin(); it != g[n_one].end(); ) {
                if (*it == n_two) {
                    edge_needed = false;
                }
                if (*it == v) {
                    it = g[n_one].erase(it);
                } else {
                    ++it;
                }
            }

            for (auto it = g[n_two].begin(); it != g[n_two].end(); ) {
                if (*it == v) {
                    it = g[n_two].erase(it);
                } else {
                    ++it;
                }
            }

            if (edge_needed) {
                g[n_one].push_back(n_two);
                g[n_two].push_back(n_one);
            }

            if (g[n_one].size() <= 2) {
                low_degree.push(n_one);
            }

            if (g[n_two].size() <= 2) {
                low_degree.push(n_two);
            }
            g[v].clear();
        }

        tr.tree.emplace_back();
        tr.nodes.emplace_back();
        int last_v = low_degree.top();
        while (g[last_v].empty()) {
            low_degree.pop();
            last_v = low_degree.top();
        }
        tr.nodes.back().insert(last_v);
        tr.nodes.back().insert(g[last_v].front());

        for (auto& node : deleted) {
            tr.tree.emplace_back();
            tr.nodes.emplace_back();
            for (int v : node) {
                tr.nodes.back().insert(v);
            }
        }

        for (int i = deleted.size() - 1; i >= 0; i--) {
            std::vector<int>& front = deleted[i];
            std::set<int> neibs;
            neibs.insert(front[1]);
            if (front.size() > 2) {
                neibs.insert(front[2]);
            }
            int node_it = 0;
            for (auto& set : tr.nodes) {
                if (node_it != i + 1 && std::includes(set.begin(), set.end(), neibs.begin(), neibs.end()) && !check_for_edge(i + 1, node_it, tr.tree)) {
                    break;
                }
                node_it++;
            }
            if (node_it < tr.nodes.size()) {
                tr.tree[node_it].push_back(i + 1);
                tr.tree[i + 1].push_back(node_it);
            }
        }
    }

    void expand_vertices() {
        std::map<std::pair<int, int>, std::vector<int>> faces;
        std::vector<std::vector<int>> vertices_in_face;
        face_getter visitor(faces, vertices_in_face);
        level_face_traversal<Graph>(embedding, visitor);

        std::vector<int> temp_outer_face = outer_face;
        for (int v : outer_face) {
            temp_outer_face.push_back(v);
        }

        for (outer_face_it = 0; outer_face_it < vertices_in_face.size(); outer_face_it++) {
            auto& face = vertices_in_face[outer_face_it];
            if (std::search(temp_outer_face.begin(), temp_outer_face.end(), face.begin(), face.end()) != temp_outer_face.end()) {
                break;
            }
        }

        std::vector<int> face_level(vertices_in_face.size(), INT16_MAX);
        face_level[outer_face_it] = 0;

        std::queue<int> Q({outer_face_it});
        while (!Q.empty()) {
            int face = Q.front();
            Q.pop();
            int prev_v = vertices_in_face[face].back();
            for (int v : vertices_in_face[face]) {
                std::vector<int>& fs = faces[std::minmax(prev_v, v)];
                int nei = fs[0] == face ? fs[1] : fs[0];
                if (face_level[nei] == INT16_MAX) {
                    Q.push(nei);
                } else {
                    face_level[face] = std::min(face_level[face], face_level[nei] + 1);
                }
                prev_v = v;
            }
        }

        std::vector<std::vector<int>> added(embedding.size());
        int size = embedding.size();
        for (int v = 0; v < size; v++) {
            if (embedding[v].size() <= 3) {
                continue;
            }

            int starting_e_it;
            int lowest_level = INT16_MAX;
            int face_num;

            Edge& prev_e = embedding[v].back();
            std::pair<int, int> prev_edge = std::minmax(prev_e.m_source, prev_e.m_target);
            for (int i = 0; i < embedding[v].size(); i++) {
                Edge& e = embedding[v][i];
                auto edge = std::minmax(e.m_source, e.m_target);
                auto &f_e = faces[edge], &f_pe = faces[prev_edge];
                auto face = std::find_if(
                    f_e.rbegin(), f_e.rend(),
                    [&] (auto f) { return std::find(f_pe.begin(), f_pe.end(), f) != f_pe.end(); });
                if (face == f_e.rend()) {
                    continue;
                }
                if (face_level[*face] < lowest_level ) {
                    lowest_level = face_level[*face];
                    face_num = *face;
                    starting_e_it = i;
                }
                prev_edge = edge;
            }

            cyclic_vector<Edge> v_emb = embedding[v];
            embedding[v].clear();
            embedding[v].push_back(v_emb[starting_e_it]);
            embedding[v].push_back(v_emb[starting_e_it + 1]);
            embedding[v].emplace_back(v, embedding.size(), nullptr);
            int last = v;

            for (int i = starting_e_it + 2; i < v_emb.size() + starting_e_it - 2; i++) {
                int curr = embedding.size();
                embedding.emplace_back();
                embedding.back().emplace_back(last, curr, nullptr);

                int nei = v_emb[i].m_source == v ? v_emb[i].m_target : v_emb[i].m_source;
                embedding.back().emplace_back(curr, nei, nullptr);
                embedding[nei][get_edge_it(v_emb[i], nei, embedding)] = Edge(curr, nei, nullptr);

                embedding.back().emplace_back(curr, curr + 1, nullptr);
                last = curr;

                expanded_vertices.push_back(v);
                added[v].push_back(curr);
            }

            int curr = embedding.size();
            embedding.emplace_back();
            embedding.back().emplace_back(last, curr, nullptr);

            int e_it = starting_e_it - 2;
            int nei = v_emb[e_it].m_source == v ? v_emb[e_it].m_target : v_emb[e_it].m_source;
            embedding.back().emplace_back(curr, nei, nullptr);
            embedding[nei][get_edge_it(v_emb[e_it], nei, embedding)] = Edge(curr, nei, nullptr);

            e_it = starting_e_it - 1;
            nei = v_emb[e_it].m_source == v ? v_emb[e_it].m_target : v_emb[e_it].m_source;
            embedding.back().emplace_back(curr, nei, nullptr);
            embedding[nei][get_edge_it(v_emb[e_it], nei, embedding)] = Edge(curr, nei, nullptr);

            expanded_vertices.push_back(v);
            added[v].push_back(curr);
        }
        faces.clear();

        Graph new_graph;
        for (int v = 0; v < embedding.size(); v++) {
            auto& v_emb = embedding[v];
            for (int i = 0; i < v_emb.size(); i++) {
                Edge edge = v_emb[i];
                int nei = edge.m_source == v ? edge.m_target : edge.m_source;
                if (v < nei) {
                    embedding[v][i] = embedding[nei][get_edge_it(edge, nei, embedding)]
                        = add_edge(edge.m_source, edge.m_target, new_graph).first;
                }
            }
        }

        graph = new_graph;
        vertices_in_face.clear();
        face_getter new_visitor(faces, vertices_in_face);
        level_face_traversal<Graph>(embedding, new_visitor);

        bool found_outer = false;
        for (outer_face_it = 0; outer_face_it < vertices_in_face.size(); outer_face_it++) {
            std::list<int> face_temp(vertices_in_face[outer_face_it].begin(), vertices_in_face[outer_face_it].end());
            for (int& v : face_temp) {
                if (expanded_vertices[v] > -1) {
                    v = expanded_vertices[v];
                }
            }

            int prev_v = face_temp.front();
            for(auto it = ++face_temp.begin(); it != face_temp.end(); ) {
                if(*it == prev_v) {
                    it = face_temp.erase(it);
                } else {
                    prev_v = *it;
                    ++it;
                }
            }

            while (face_temp.front() == face_temp.back()) {
                face_temp.pop_back();
            }

            if (std::search(temp_outer_face.begin(), temp_outer_face.end(), face_temp.begin(), face_temp.end())
                != temp_outer_face.end()) {
                break;
            }
        }

        outer_face.clear();
        for (int v : vertices_in_face[outer_face_it]) {
            outer_face.push_back(v);
        }
    }

    void create_spanning_forest(std::vector<std::vector<NetworKit::Edge>> &outer_edges) {
        std::vector<std::vector<int>> vertices_on_level(k + 1);

        for (int v = 0; v < vertex_level.size(); v++) {
            vertices_on_level[vertex_level[v]].push_back(v);
        }

        for (int level = 1; level <= k; level++) {
            std::vector<cyclic_vector<int>>& level_graph = level_graphs[level];
            std::vector<NetworKit::Edge>& edges = outer_edges[level];
            for (const auto &[v, w] : edges) {
                if (!check_for_edge(v, w, level_graph)) {
                    level_graph[v].push_back(w), level_graph[w].push_back(v);
                }
            }
            for (int v : vertices_on_level[level]) {
                for (auto& edge : embedding[v]) {
                    int w = edge.m_source == v ? edge.m_target : edge.m_source;
                    if (vertex_level[w] == level - 1) {
                        level_graph[v].push_back(w), level_graph[w].push_back(v);
                        edges.emplace_back(v, w);
                    }
                }
            }
            std::vector<int> visited(embedding.size(), 0);
            cycle_num = 0;
            for (int v : vertices_on_level[level]) {
                if (visited[v] == 0) {
                    cycle_dfs(v, -1, level, visited);
                }
            }
        }

        std::vector<bool> vertex_in_tree(embedding.size());
        for (int v : vertices_on_level.back()) {
            for (auto& edge : embedding[v]) {
                int w = edge.m_source == v ? edge.m_target : edge.m_source;
                if (vertex_level[v] == k && vertex_level[w] == k && !check_for_edge(v, w, level_graphs[k])) {
                    spanning_forest[v].push_back(w);
                    vertex_in_tree[v] = true;
                }
            }
        }

        // 0 = gamma, 1 = alpha, 2 = beta
        std::vector<int> category(embedding.size());
        for (int level = k; level > 0; level--) {
            for (int v : vertices_on_level[k]) {
                category[v] = 0;
            }
            auto &cycles = level_cycles[level], &level_graph = level_graphs[level];
            for (auto& cycle : cycles) {
                for (int v : cycle) {
                    if (level_graph[v].size() == 3) {
                        category[v] = 1;
                    } else if (vertex_in_tree[v]) {
                        category[v] = 2;
                    }
                }
            }
            for (auto &[v, w] : outer_edges[level]) {
                if (category[v] == 0 && category[w] == 0 && !check_for_edge(v, w, spanning_forest)) {
                    spanning_forest[v].push_back(w), spanning_forest[w].push_back(v);
                    vertex_in_tree[v] = vertex_in_tree[w] = true;
                }
            }

            for (auto& cycle : cycles) {
                std::vector<std::pair<int, int>> edges_to_add;

                int starting_alpha_it = -1;
                for (int i = 0; i < cycle.size(); i++) {
                    int v = cycle[i];

                    if (category[v] == 1) {
                        for (int w : level_graph[v]) {
                            if (w != cycle[i - 1] && w != cycle[i + 1]) {
                                int x = cycle[-1];
                                edges_to_add.emplace_back(v, w);
                                vertex_in_tree[v] = vertex_in_tree[w] = true;
                            }
                        }
                        starting_alpha_it = i;
                    }
                }

                if (starting_alpha_it == -1) {
                    for (int i = 0; i < cycle.size(); i++) {
                        int v = cycle[i];

                        if (category[v] == 0) {
                            starting_alpha_it = i;
                        }
                    }
                }

                std::set<int> leaves;
                if (starting_alpha_it == -1) {
                    starting_alpha_it = 0;
                    get_leaves(cycle[0], -1, level, leaves);
                }

                int prev_v = cycle[starting_alpha_it];
                for (int v_it = starting_alpha_it + 1; v_it < starting_alpha_it + cycle.size(); v_it++) {
                    int v = cycle[v_it];

                    if (category[v] == 2 && leaves.find(v) != leaves.end()) {
                        prev_v = v;
                        continue;
                    }

                    if (category[v] == 2) {
                        get_leaves(v, -1, level, leaves);
                    }

                    edges_to_add.emplace_back(v, prev_v);
                    vertex_in_tree[v] = vertex_in_tree[prev_v] = true;
                    prev_v = v;
                }

                for (auto &[u, v] : edges_to_add) {
                    if (!check_for_edge(u, v, spanning_forest)) {
                        spanning_forest[u].push_back(v), spanning_forest[v].push_back(u);
                    }
                }

                if (check_for_edge(cycle[starting_alpha_it], cycle[starting_alpha_it - 1], spanning_forest)) {
                    int v = cycle[starting_alpha_it], w = cycle[starting_alpha_it - 1];
                    spanning_forest[v].erase(
                        std::remove(spanning_forest[v].begin(), spanning_forest[v].end(), w),
                        spanning_forest[v].end());
                    spanning_forest[w].erase(
                        std::remove(spanning_forest[w].begin(), spanning_forest[w].end(), v),
                        spanning_forest[w].end());
                }
            }
        }
    }

    struct node{
        std::pair<int, Edge> p;
        std::vector<std::pair<int, Edge>> children;
        std::vector<Edge> face;
    };

    void tr_dfs(int v, int p, std::vector<Edge>& edges, std::vector<node>& tree) {
        Edge missing_e;
        for (auto c : tree[v].children) {
            if (c.first == p) {
                missing_e = c.second;
            } else {
                std::vector<Edge> temp_edges;
                tr_dfs(c.first, v, temp_edges, tree);
                for (auto& e : temp_edges) {
                    edges.push_back(e);
                }
            }
        }
        if (p == -1) {
            return;
        }
        for (auto& e : tree[v].face) {
            edges.push_back(e);
        }
        int end = missing_e.m_source, other_end = missing_e.m_target;
        for (auto &e : edges) {
            tr[e].insert(end);
            if (e.m_source != other_end) {
                tr[e.m_source].insert(end);
            }
            if (e.m_target != other_end) {
                tr[e.m_target].insert(end);
            }
        }
    }

    void build_tree_decomposition() {
        for (int v = 0; v < embedding.size(); v++) {
            tr.nodes.emplace_back();
            tr.tree.emplace_back();
            tr[v].insert(v);
        }
        std::map<std::pair<int, int>, std::vector<int>> faces;
        std::vector<std::vector<int>> vertices_in_face;
        face_getter visitor(faces, vertices_in_face);
        level_face_traversal<Graph>(embedding, visitor);
        std::vector<node> face_tree(vertices_in_face.size());
        for (int face = 0; face < vertices_in_face.size(); face++) {
            int prev_v = vertices_in_face[face].back();
            for (int v : vertices_in_face[face]) {
                if (check_for_edge(prev_v, v, spanning_forest)) {
                    face_tree[face].face.push_back(Edge(prev_v, v, nullptr));
                }
                prev_v = v;
            }
        }
        boost::graph_traits<Graph>::edge_iterator ei, ei_end;
        for(boost::tie(ei, ei_end) = edges(graph); ei != ei_end; ++ei) {
            int v = std::min(ei->m_source, ei->m_target);
            int w = std::max(ei->m_source, ei->m_target);
            std::vector<int>& fs = faces[std::make_pair(v, w)];
            if (!check_for_edge(v, w, spanning_forest)) {
                face_tree[fs[0]].children.push_back(std::make_pair(fs[1], *ei));
                face_tree[fs[1]].children.push_back(std::make_pair(fs[0], *ei));
            } else {
                int s = tr.nodes.size();
                tr.nodes.emplace_back();
                tr.tree.emplace_back();
                tr[s].insert(v), tr[s].insert(w);
                std::pair<int, int> e = std::minmax(ei->m_source, ei->m_target);
                tr.edge_it[e] = s;
                tr.tree[s].push_back(v), tr.tree[v].push_back(s);
                tr.tree[s].push_back(w), tr.tree[w].push_back(s);
            }
        }
        std::vector<Edge> edges;
        tr_dfs(outer_face_it, -1, edges, face_tree);
    }

    void contract_vertices(){
        for (auto& node : tr.nodes) {
            std::vector<int> insert_list;
            for(auto it = node.begin(); it != node.end(); ) {
                if(expanded_vertices[*it] != -1) {
                    insert_list.push_back(expanded_vertices[*it]);
                    it = node.erase(it);
                } else {
                    ++it;
                }
            }

            for (int v : insert_list) {
                node.insert(v);
            }
        }
    }

public:
    create_tree_decomposition(Graph& g, std::vector<cyclic_vector<Edge>>& emb, std::vector<int> out_face)
            : graph(g), embedding(emb), outer_face(out_face), expanded_vertices(embedding.size(), -1), vertex_level(embedding.size()) {
        auto [k_prim, _] = name_levels(embedding, out_face, vertex_level);
        this->k = k_prim;
//        if (k == 1) {
//            outerplanar();
//            return;
//        }
        expand_vertices();
        vertex_level.resize(embedding.size());
        spanning_forest.resize(embedding.size());
        auto [k_bis, outer_edges] = name_levels(embedding, outer_face, vertex_level);
        this->k = k_bis;
        outer_edges.emplace_back();
        level_graphs = std::vector<std::vector<cyclic_vector<int>>>(k + 1, std::vector<cyclic_vector<int>>(embedding.size()));
        level_cycles.resize(k + 1);
        create_spanning_forest(outer_edges);
        build_tree_decomposition();
        contract_vertices();
    }

    tree_decomposition get_tree_decomposition() {
        return tr;
    }
};

void get_tree_decomposition(Graph& g, PlanarEmbedding& emb, std::vector<int> outer_face, tree_decomposition& tr) {
    create_tree_decomposition ctr(g, emb, outer_face);
    tr = ctr.get_tree_decomposition();
}
