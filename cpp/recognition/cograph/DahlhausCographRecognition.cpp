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
        is_cograph = check_cotree(T);
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
        std::vector<bool> &is_neighbour) {
        for (auto u : G.neighborRange(v)) {
            if (is_neighbour[u] || component[u] != -1) {
                continue;
            }
            component[u] = component[v];
            dfs(u, G, component, is_neighbour);
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
        if (vertex_type == 0) {
            u2->AddChild(T2.root);
        } else {
            auto u1 = T.Add(Type::ZERO_ONE, 0);
            u2->AddChild(u1);
            u1->AddChild(T2.root);
        }
    }

    CoTree &DahlhausCographRecognition::build_cotree(NetworKit::Graph G,
        std::vector<int> real_number_of_node) {  // should return cotree reference
        int n = G.numberOfNodes();
        auto v = *G.nodeRange().begin();
        CoTree &T = save.emplace_back();
        T.ReserveSpace(3 * n);
        auto V = T.Add(Type::VERTEX, real_number_of_node[v]);
        pointer[real_number_of_node[v]] = V;
        T.root = V;
        if (n == 1) {
            return T;
        }
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
        int component_number = 0;
        for (auto u : not_neighbours) {
            if (component[u] == -1) {
                component[u] = component_number++;
                dfs(u, G, component, is_neighbour);
            }
        }
        std::vector<std::vector<int>> components(component_number), gamma(n);
        for (int i = 0; i < n; i++) {
            if (is_neighbour[i] || i == v) {
                continue;
            }
            components[component[i]].push_back(i);
        }
        for (int i = 0; i < n; i++) {
            if (is_neighbour[i] || i == v) {
                continue;
            }
            for (auto u : G.neighborRange(i)) {
                if (!is_neighbour[u] && component[u] == component[i]) {
                    continue;
                }
                gamma[i].push_back(u);
            }
        }
        std::sort(components.begin(), components.end(),
            [&gamma](std::vector<int> &a, std::vector<int> &b)
            {return gamma[a[0]].size() > gamma[b[0]].size();
        });
        for (int i = 0; i < component_number; i++) {
            for (int j = 0; j < components[i].size(); j++) {
                component[components[i][j]] = i;
            }
        }
        std::vector<std::vector<int>> gamma_difference(component_number + 1);
        std::vector<int> last_position_where_met(n, -1);
        for (int i = 0; i < n; i++) {
            if (is_neighbour[i] || i == v) {
                continue;
            }
            for (auto a : gamma[i]) {
                last_position_where_met[a] = std::max(last_position_where_met[a], component[i]);
            }
        }
        for (int i = 0; i < n; i++) {
            if (!is_neighbour[i] || i == v) {
                continue;
            }
            gamma_difference[1 + last_position_where_met[i]].push_back(i);
        }
        std::vector<int> fake_number_of_node(n, -1);
        for (int i = 0; i <= component_number; i++) {
            add(1, T, gamma_difference[i], fake_number_of_node, G, real_number_of_node);
            if (i == component_number) {
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
