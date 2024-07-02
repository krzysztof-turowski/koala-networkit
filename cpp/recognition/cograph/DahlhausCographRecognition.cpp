/*
 * DahlhausCographRecognition.cpp
 *
 *  Created on: 2024
 *      Author: fixikmila
 */
#include <graph/GraphTools.hpp>

#include "recognition/CographRecognition.hpp"
#include "recognition/CoTree.hpp"

namespace Koala {

    int current_time = 0;
    std::vector<CoNode *> dfs_list;
    int number_of_edges_according_to_cotree = 0;
    int maximum_depth = 0;

    void dfsCoTree(CoNode *v, int depth = 0) {
        maximum_depth = std::max(maximum_depth, depth);
        dfs_list.push_back(v);
        v->time_in = current_time++;
        if (v->type == Type::VERTEX) {
            v->number_of_vertices_in_subtree = 1;
            v->time_out = current_time++;
            return;
        }
        CoNode *child = v->first_child;
        int sum = 0;
        while (child != nullptr) {
            dfsCoTree(child, depth + 1);
            if (v->number == 1) {
                number_of_edges_according_to_cotree += child->number_of_vertices_in_subtree * sum;
            }
            sum += child->number_of_vertices_in_subtree;
            v->number_of_vertices_in_subtree += child->number_of_vertices_in_subtree;
            child = child->next;
        }
        v->time_out = current_time++;
    }

    bool descendant(CoNode *u, CoNode *v) {
        return u->time_in >= v->time_in && u->time_out <= v->time_out;
    }

    CoNode *lca(CoNode *u, CoNode *v, int logarithm) {
        if (descendant(u, v)) {
            return v;
        }
        if (descendant(v, u)) {
            return u;
        }
        for (int i = logarithm - 1; i >= 0; i--) {
            if (!descendant(v, u->get_up[i])) {
                u = u->get_up[i];
            }
        }
        return u->get_up[0];
    }

    void DahlhausCographRecognition::run() {
        is_cograph = true;
        hasRun = true;
        std::vector<NetworKit::node> nodes;
        for (auto u : graph.nodeRange()) {
            nodes.push_back(u);
        }
        std::vector<int> real_number_of_node(graph.numberOfNodes());
        for (int i = 0; i < graph.numberOfNodes(); i++) {
            real_number_of_node[i] = i;
        }
        // reserve space for cotree and do in the same way as conodes are stored
        save.reserve(nodes.size() * 4);
        pointer.resize(nodes.size());
        auto &T = build_cotree(graph, real_number_of_node);
        if (is_cograph) {
            is_cograph = check_cotree(T);
        }
        for (auto &t : save) {
            t.Clear();
        }
        save.clear();
    }

    bool DahlhausCographRecognition::isCograph() const {
        assureFinished();
        return is_cograph;
    }

    void dfs(NetworKit::node v, NetworKit::Graph &G, std::vector<int> &component,
        std::vector<bool> &is_in_vec) {
        for (auto u : G.neighborRange(v)) {
            if (!is_in_vec[u] || component[u] != -1) {
                continue;
            }
            component[u] = component[v];
            dfs(u, G, component, is_in_vec);
        }
    }

    inline NetworKit::Graph
    build_graph(std::vector<int> &nodes, NetworKit::Graph &G,
        std::vector<int> &fake_number_of_node) {
        NetworKit::Graph h;
        for (auto u : nodes) {
            h.addNode();
        }
        for (auto u : nodes) {
            for (auto v : G.neighborRange(u)) {
                if (fake_number_of_node[v] != -1 && u < v) {
                    h.addEdge(fake_number_of_node[u], fake_number_of_node[v]);
                }
            }
        }
        return h;
    }

    inline void DahlhausCographRecognition::add(int vertex_type, CoTree &T, std::vector<int> &vec,
        std::vector<int> &fake_number_of_node, NetworKit::Graph &G,
        std::vector<int> &real_number_of_node) {
        auto u2 = T.Add(Type::ZERO_ONE, vertex_type);
        u2->AddChild(T.root);
        T.root = u2;
        if (vec.empty()) {
            return;
        }
        for (int j = 0; j < vec.size(); j++) {
            fake_number_of_node[vec[j]] = j;
        }
        auto C = build_graph(vec, G, fake_number_of_node);
        for (int j = 0; j < vec.size(); j++) {
            fake_number_of_node[vec[j]] = -1;
        }
        std::vector<int> new_real_number_of_node(vec.size());
        for (int j = 0; j < vec.size(); j++) {
            new_real_number_of_node[j] = real_number_of_node[vec[j]];
        }
        CoTree &T2 = build_cotree(C, new_real_number_of_node);
        if (!is_cograph) {
            return;
        }
        if (vertex_type == 0) {
            u2->AddChild(T2.root);
        } else {
            auto u1 = T.Add(Type::ZERO_ONE, 0);
            u2->AddChild(u1);
            u1->AddChild(T2.root);
        }
    }

    std::vector<std::vector<int>>
    compute_connected_components(std::vector<int> &vec, std::vector<int> &component,
        std::vector<bool> &is_in_vec, NetworKit::Graph &G) {
        int component_number = 0;
        for (auto u : vec) {
            if (component[u] == -1) {
                component[u] = component_number++;
                dfs(u, G, component, is_in_vec);
            }
        }
        std::vector<std::vector<int>> components(component_number);
        for (int i = 0; i < is_in_vec.size(); i++) {
            if (!is_in_vec[i]) {
                continue;
            }
            components[component[i]].push_back(i);
        }
        return components;
    }

    std::vector<std::vector<int>>
    compute_gamma(std::vector<bool> &is_in_vec, NetworKit::Graph &G, std::vector<int> &component) {
        std::vector<std::vector<int>> gamma(is_in_vec.size());
        for (int i = 0; i < is_in_vec.size(); i++) {
            if (!is_in_vec[i]) {
                continue;
            }
            for (auto u : G.neighborRange(i)) {
                if (is_in_vec[u] && component[u] == component[i]) {
                    continue;
                }
                gamma[i].push_back(u);
            }
        }
        return gamma;
    }

    std::vector<std::vector<int>>
    compute_components_sorted(int n, std::vector<std::vector<int>> &components,
        std::vector<std::vector<int>> &gamma) {
        std::vector<std::vector<int>> count_sort(n), components_sorted;
        for (int i = 0; i < components.size(); i++) {
            count_sort[gamma[components[i][0]].size()].push_back(i);
        }
        for (int i = n - 1; i >= 0; i--) {
            for (auto value : count_sort[i]) {
                components_sorted.push_back(components[value]);
            }
        }
        return components_sorted;
    }

    void recompute_component(std::vector<std::vector<int>> &components,
                             std::vector<int> &component) {
        for (int i = 0; i < components.size(); i++) {
            for (int j = 0; j < components[i].size(); j++) {
                component[components[i][j]] = i;
            }
        }
    }

    std::vector<std::vector<int>>
    compute_gamma_difference(std::vector<std::vector<int>> &components, std::vector<int> &component,
        std::vector<std::vector<int>> &gamma, std::vector<bool> &is_in_vec,
        std::vector<bool> &is_in_new_vec) {
        std::vector<std::vector<int>> gamma_difference(components.size() + 1);
        int n = is_in_vec.size();
        std::vector<int> last_position_where_met(n, -1);
        for (int i = 0; i < n; i++) {
            if (!is_in_vec[i]) {
                continue;
            }
            for (auto a : gamma[i]) {
                last_position_where_met[a] = std::max(last_position_where_met[a], component[i]);
            }
        }
        for (int i = 0; i < n; i++) {
            if (is_in_new_vec[i]) {
                continue;
            }
            gamma_difference[1 + last_position_where_met[i]].push_back(i);
        }
        return gamma_difference;
    }

    void reverse_cotree(CoTree &T, CoNode *v) {
        if (v->type == Type::ZERO_ONE) {
            v->number ^= 1;
        }
        auto u = v->first_child;
        while (u != nullptr) {
            reverse_cotree(T, u);
            u = u->next;
        }
    }

    void DahlhausCographRecognition::big_component(CoTree &T, NetworKit::Graph &G,
        std::vector<int> &vec, std::vector<int> &real_number_of_node) {
        NetworKit::Graph GC = Koala::GraphTools::toComplement(G);
        int n = GC.numberOfNodes();
        std::vector<int> component(n, -1), fake_number_of_node(n, -1);
        std::vector<bool> is_in_vec(n);
        for (auto u : vec) {
            is_in_vec[u] = true;
        }
        auto components = compute_connected_components(vec, component, is_in_vec, GC);
        for (auto c : components) {
            if (c.size() * A > 2 * n + A) {
                is_cograph = false;
                return;
            }
            for (int j = 0; j < c.size(); j++) {
                fake_number_of_node[c[j]] = j;
            }
            auto GI = build_graph(c, GC, fake_number_of_node);  // induced subgraph by c
            for (int j = 0; j < c.size(); j++) {
                fake_number_of_node[c[j]] = -1;
            }
            std::vector<int> new_real_number_of_node(c.size());
            for (int j = 0; j < c.size(); j++) {
                new_real_number_of_node[j] = real_number_of_node[c[j]];
            }
            auto TI = build_cotree(GI, new_real_number_of_node);
            if (!is_cograph) {
                return;
            }
            reverse_cotree(TI, TI.root);  // reverse TI
            T.root->AddChild(TI.root);
        }
    }

    void
    DahlhausCographRecognition::high_low_case(CoTree &T, NetworKit::Graph &G,
        std::vector<int> &real_number_of_node) {
        if (!is_cograph) {
            return;
        }
        int n = G.numberOfNodes();
        std::vector<int> degree(n);
        for (auto u : G.nodeRange()) {
            for (auto q : G.neighborRange(u)) {
                degree[u]++;
            }
        }
        auto V = T.Add(Type::ZERO_ONE, 0);
        T.root = V;
        // compute low components
        // sort and compute gamma difference
        // 0-low components 1-high components
        // if high component is big, then call big_components
        std::vector<bool> is_in_vec(n);
        std::vector<int> vec, component(n, -1);
        for (auto u : G.nodeRange()) {
            if (degree[u] * A <= n) {
                is_in_vec[u] = true;
                vec.push_back(u);
            }
        }
        auto components = compute_connected_components(vec, component, is_in_vec, G);
        auto gamma = compute_gamma(is_in_vec, G, component);
        components = compute_components_sorted(n, components, gamma);
        recompute_component(components, component);
        auto gamma_difference =
                compute_gamma_difference(components, component, gamma, is_in_vec, is_in_vec);
        std::vector<int> fake_number_of_node(n, -1);

        for (int i = 0; i <= components.size(); i++) {
            bool special_case_big_component = true;
            if (gamma_difference[i].size() * A > (A - 1) * n) {
                std::vector<bool> is_in_gamma_difference(n);
                for (auto u : gamma_difference[i]) {
                    is_in_gamma_difference[u] = true;
                }
                for (auto u : gamma_difference[i]) {
                    int sum = 0;
                    for (auto v : G.neighborRange(u)) {
                        if (!is_in_gamma_difference[v]) {
                            continue;
                        }
                        sum++;
                    }
                    sum = int(gamma_difference[i].size()) - 1 - sum;
                    if (sum * A >= n) {
                        special_case_big_component = false;
                        break;
                    }
                }
            } else {
                special_case_big_component = false;
            }
            if (special_case_big_component) {
                std::vector<int> empty;
                add(1, T, empty, fake_number_of_node, G, real_number_of_node);
                big_component(T, G, gamma_difference[i], real_number_of_node);
                if (!is_cograph) {
                    return;
                }
            } else {
                add(1, T, gamma_difference[i], fake_number_of_node, G, real_number_of_node);
            }
            if (i == components.size()) {
                break;
            }
            if (components[i].size() * A > 2 * n + A) {
                is_cograph = false;
                return;
            }
            add(0, T, components[i], fake_number_of_node, G, real_number_of_node);
        }
    }

    CoTree &DahlhausCographRecognition::build_cotree(NetworKit::Graph G,
        std::vector<int> real_number_of_node) {  // should return cotree reference
        int n = G.numberOfNodes();
        CoTree &T = save.emplace_back();
        T.ReserveSpace(3 * n);
        if (n == 1) {
            auto v = *G.nodeRange().begin();
            auto V = T.Add(Type::VERTEX, real_number_of_node[v]);
            pointer[real_number_of_node[v]] = V;
            T.root = V;
            return T;
        }
        if (!is_cograph) {
            return T;
        }

        int v = -1;
        for (auto u : G.nodeRange()) {
            int size = 0;
            for (auto q : G.neighborRange(u)) {
                size++;
            }
            if (A * size >= n && size * A <= (A - 1) * n) {
                v = u;
                break;
            }
        }
        if (v == -1) {
            high_low_case(T, G, real_number_of_node);
            return T;
        }
        auto V = T.Add(Type::VERTEX, real_number_of_node[v]);
        pointer[real_number_of_node[v]] = V;
        T.root = V;
        std::vector<bool> is_neighbour(n);
        std::vector<int> not_neighbours;
        for (auto u : G.neighborRange(v)) {
            is_neighbour[u] = true;
        }
        for (unsigned int i = 0; i < n; i++) {
            if (i == v) {
                continue;
            }
            if (!is_neighbour[i]) {
                not_neighbours.push_back(i);
            }
        }
        std::vector<int> component(n, -1);
        std::vector<bool> is_in_vec(n, true);
        for (int i = 0; i < n; i++) {
            if (is_neighbour[i] || i == v) {
                is_in_vec[i] = false;
            }
        }
        std::vector<std::vector<int>> components =
                compute_connected_components(not_neighbours, component, is_in_vec, G),
                gamma = compute_gamma(is_in_vec, G, component);
        components = compute_components_sorted(n, components, gamma);
        recompute_component(components, component);
        std::vector<bool> is_in_new_vec = is_in_vec;
        is_in_new_vec[v] = true;
        auto gamma_difference =
                compute_gamma_difference(components, component, gamma, is_in_vec, is_in_new_vec);
        std::vector<int> fake_number_of_node(n, -1);
        for (int i = 0; i <= components.size(); i++) {
            add(1, T, gamma_difference[i], fake_number_of_node, G, real_number_of_node);
            if (i == components.size()) {
                break;
            }
            add(0, T, components[i], fake_number_of_node, G, real_number_of_node);
        }
        return T;
    }

    bool DahlhausCographRecognition::check_cotree(CoTree T) {
        current_time = 0;
        dfs_list.clear();
        number_of_edges_according_to_cotree = 0;
        maximum_depth = 0;
        dfsCoTree(T.root);
        for (auto u : dfs_list) {
            if (u->parent == nullptr) {
                u->get_up[0] = u;
            } else {
                u->get_up[0] = u->parent;
            }
        }
        int i = 1;
        for (int number = 2; number <= maximum_depth; i++, number <<= 1) {
            for (auto u : dfs_list) {
                u->get_up[i] = u->get_up[i - 1]->get_up[i - 1];
            }
        }
        for (auto [u, v] : graph.edgeRange()) {
            auto l = lca(pointer[u], pointer[v], i);
            if (l->type != Type::ZERO_ONE || l->number != 1) {
                return false;
            }
        }
        return number_of_edges_according_to_cotree == graph.numberOfEdges();
    }

}  // namespace Koala
