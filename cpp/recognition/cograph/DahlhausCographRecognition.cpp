/*
 * DahlhausCographRecognition.cpp
 *
 *  Created on: 07.06.2024
 *      Author: fixikmila
 */

#include <algorithm>
#include <numeric>
#include <vector>

#include <graph/GraphTools.hpp>
#include <recognition/CographRecognition.hpp>
#include <structures/Cotree.hpp>

namespace Koala {

void DahlhausCographRecognition::run() {
    is_cograph = State::COGRAPH;
    hasRun = true;
    std::vector<NetworKit::node> nodes;
    for (auto u : graph.nodeRange()) {
        nodes.push_back(u);
    }
    std::vector<int> real_index(graph.numberOfNodes());
    std::iota(real_index.begin(), real_index.end(), 0);
    cotree.clear();
    cotree.reserve(nodes.size() * 4);
    pointer.assign(nodes.size(), NetworKit::none);
    auto root = build_cotree(cotree, graph, real_index);
    cotree.setRoot(root);
    if (is_cograph == State::COGRAPH) {
        is_cograph = check_cotree(cotree) ? State::COGRAPH : State::NOT_COGRAPH;
    }
    cotree.clear();
}

bool descendant(Cotree &T, NetworKit::count u, NetworKit::count v) {
    return T.getNode(u).time_in >= T.getNode(v).time_in
        && T.getNode(u).time_out <= T.getNode(v).time_out;
}

NetworKit::count lca(Cotree &T, NetworKit::count u, NetworKit::count v, int logarithm) {
    if (descendant(T, u, v)) {
        return v;
    }
    if (descendant(T, v, u)) {
        return u;
    }
    for (int i = logarithm - 1; i >= 0; i--) {
        if (!descendant(T, v, T.getNode(u).get_up[i])) {
            u = T.getNode(u).get_up[i];
        }
    }
    return T.getNode(u).get_up[0];
}

void dfs(
        NetworKit::node v, NetworKit::Graph &G, std::vector<int> &component,
        std::vector<bool> &is_in_vec) {
    for (auto u : G.neighborRange(v)) {
        if (!is_in_vec[u] || component[u] != -1) {
            continue;
        }
        component[u] = component[v];
        dfs(u, G, component, is_in_vec);
    }
}

inline NetworKit::Graph build_graph(
        std::vector<int> &nodes, NetworKit::Graph &G, std::vector<int> &fake_index) {
    NetworKit::Graph h(nodes.size());
    for (auto u : nodes) {
        for (auto v : G.neighborRange(u)) {
            if (fake_index[v] != -1 && static_cast<NetworKit::node>(u) < v) {
                h.addEdge(fake_index[u], fake_index[v]);
            }
        }
    }
    return h;
}

inline void DahlhausCographRecognition::add(
        int vertex_type, Cotree &T, std::vector<int> &vec,
        std::vector<int> &fake_index, NetworKit::Graph &G, std::vector<int> &real_index) {
    auto u2 = T.add(
        vertex_type == 0 ? NodeType::UNION_NODE : NodeType::COMPLEMENT_NODE, vertex_type);
    T.addChild(u2, T.getRoot());
    T.setRoot(u2);
    if (vec.empty()) {
        return;
    }
    int fake = 0;
    for (auto v : vec) {
        fake_index[v] = fake++;
    }
    auto C = build_graph(vec, G, fake_index);
    for (auto v : vec) {
        fake_index[v] = -1;
    }
    std::vector<int> new_real_index(vec.size());
    for (std::size_t j = 0; j < vec.size(); j++) {
        new_real_index[j] = real_index[vec[j]];
    }
    auto root = T.getRoot();
    auto subtree_root = build_cotree(T, C, new_real_index);
    if (is_cograph != State::COGRAPH) {
        return;
    }
    T.setRoot(root);
    if (vertex_type == 0) {
        T.addChild(u2, subtree_root);
    } else {
        auto u1 = T.add(NodeType::UNION_NODE, 0);
        T.addChild(u2, u1);
        T.addChild(u1, subtree_root);
    }
    T.setRoot(u2);
}

std::vector<std::vector<int>> compute_connected_components(
        std::vector<int> &vec, std::vector<int> &component, std::vector<bool> &is_in_vec,
        NetworKit::Graph &G) {
    int component_number = 0;
    for (auto u : vec) {
        if (component[u] == -1) {
            component[u] = component_number++;
            dfs(u, G, component, is_in_vec);
        }
    }
    std::vector<std::vector<int>> components(component_number);
    for (std::size_t i = 0; i < is_in_vec.size(); i++) {
        if (!is_in_vec[i]) {
            continue;
        }
        components[component[i]].push_back(i);
    }
    return components;
}

std::vector<std::vector<int>> compute_gamma(
        std::vector<bool> &is_in_vec, NetworKit::Graph &G, std::vector<int> &component) {
    std::vector<std::vector<int>> gamma(is_in_vec.size());
    for (std::size_t i = 0; i < is_in_vec.size(); i++) {
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

std::vector<std::vector<int>> compute_components_sorted(
        NetworKit::count n, std::vector<std::vector<int>> &components,
        std::vector<std::vector<int>> &gamma) {
    std::vector<std::vector<int>> count_sort(n), components_sorted;
    for (std::size_t i = 0; i < components.size(); i++) {
        count_sort[gamma[components[i][0]].size()].push_back(i);
    }
    for (int i = n - 1; i >= 0; i--) {
        for (auto value : count_sort[i]) {
            components_sorted.push_back(components[value]);
        }
    }
    return components_sorted;
}

void recompute_component(
        std::vector<std::vector<int>> &components, std::vector<int> &component) {
    for (std::size_t i = 0; i < components.size(); i++) {
        for (std::size_t j = 0; j < components[i].size(); j++) {
            component[components[i][j]] = i;
        }
    }
}

std::vector<std::vector<int>>
compute_gamma_difference(
        std::vector<std::vector<int>> &components, std::vector<int> &component,
        std::vector<std::vector<int>> &gamma, std::vector<bool> &is_in_vec,
        std::vector<bool> &is_in_new_vec) {
    std::vector<std::vector<int>> gamma_difference(components.size() + 1);
    std::size_t n = is_in_vec.size();
    std::vector<int> last_position_where_met(n, -1);
    for (std::size_t i = 0; i < n; i++) {
        if (!is_in_vec[i]) {
            continue;
        }
        for (auto a : gamma[i]) {
            last_position_where_met[a] = std::max(last_position_where_met[a], component[i]);
        }
    }
    for (std::size_t i = 0; i < n; i++) {
        if (is_in_new_vec[i]) {
            continue;
        }
        gamma_difference[1 + last_position_where_met[i]].push_back(i);
    }
    return gamma_difference;
}

void reverse_cotree(Cotree &T, NetworKit::count v) {
    if (T.getNode(v).type != NodeType::LEAF) {
        T.getNode(v).number ^= 1;
        T.getNode(v).type = T.getNode(v).number == 0
            ? NodeType::UNION_NODE : NodeType::COMPLEMENT_NODE;
    }
    auto u = T.getNode(v).first_child;
    while (u != NetworKit::none) {
        reverse_cotree(T, u);
        u = T.getNode(u).next_sibling;
    }
}

void DahlhausCographRecognition::big_component(
        Cotree &T, NetworKit::Graph &G, std::vector<int> &vec, std::vector<int> &real_index) {
    NetworKit::Graph GC = Koala::GraphTools::toComplement(G);
    NetworKit::count n = GC.numberOfNodes();
    std::vector<int> component(n, -1), fake_index(n, -1);
    std::vector<bool> is_in_vec(n);
    for (auto u : vec) {
        is_in_vec[u] = true;
    }
    auto components = compute_connected_components(vec, component, is_in_vec, GC);
    for (auto c : components) {
        if (c.size() * A > 2 * n + A) {
            is_cograph = State::NOT_COGRAPH;
            return;
        }
        for (std::size_t j = 0; j < c.size(); j++) {
            fake_index[c[j]] = j;
        }
        auto GI = build_graph(c, GC, fake_index);  // induced subgraph by c
        for (std::size_t j = 0; j < c.size(); j++) {
            fake_index[c[j]] = -1;
        }
        std::vector<int> new_real_index(c.size());
        for (std::size_t j = 0; j < c.size(); j++) {
            new_real_index[j] = real_index[c[j]];
        }
        auto root = T.getRoot();
        auto subtree_root = build_cotree(T, GI, new_real_index);
        if (is_cograph != State::COGRAPH) {
            return;
        }
        T.setRoot(root);
        reverse_cotree(T, subtree_root);  // reverse TI
        T.addChild(T.getRoot(), subtree_root);
    }
}

void DahlhausCographRecognition::high_low_case(
        Cotree &T, NetworKit::Graph &G, std::vector<int> &real_index) {
    if (is_cograph != State::COGRAPH) {
        return;
    }
    NetworKit::count n = G.numberOfNodes();
    std::vector<NetworKit::count> degree(n);
    for (auto u : G.nodeRange()) {
        degree[u] = G.degree(u);
    }
    auto V = T.add(NodeType::UNION_NODE, 0);
    T.setRoot(V);
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
    std::vector<int> fake_index(n, -1);

    for (std::size_t i = 0; i <= components.size(); i++) {
        bool special_case_big_component = true;
        if (gamma_difference[i].size() * A > (A - 1) * n) {
            std::vector<bool> is_in_gamma_difference(n);
            for (auto u : gamma_difference[i]) {
                is_in_gamma_difference[u] = true;
            }
            for (auto u : gamma_difference[i]) {
                NetworKit::count sum = 0;
                for (auto v : G.neighborRange(u)) {
                    if (!is_in_gamma_difference[v]) {
                        continue;
                    }
                    sum++;
                }
                sum = gamma_difference[i].size() - 1 - sum;
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
            add(1, T, empty, fake_index, G, real_index);
            big_component(T, G, gamma_difference[i], real_index);
            if (is_cograph != State::COGRAPH) {
                return;
            }
        } else {
            add(1, T, gamma_difference[i], fake_index, G, real_index);
        }
        if (i == components.size()) {
            break;
        }
        if (components[i].size() * A > 2 * n + A) {
            is_cograph = State::NOT_COGRAPH;
            return;
        }
        add(0, T, components[i], fake_index, G, real_index);
    }
}

NetworKit::count DahlhausCographRecognition::build_cotree(
        Cotree &T,
        NetworKit::Graph G, std::vector<int> real_index) {  // should return cotree reference
    NetworKit::count n = G.numberOfNodes();
    T.reserve(3 * n);
    if (n == 1) {
        auto v = *G.nodeRange().begin();
        auto V = T.add(NodeType::LEAF, real_index[v]);
        pointer[real_index[v]] = V;
        T.setRoot(V);
        return V;
    }
    if (is_cograph != State::COGRAPH) {
        return T.getRoot();
    }

    NetworKit::node v = NetworKit::none;
    for (auto u : G.nodeRange()) {
        auto size = G.degree(u);
        if (A * size >= n && size * A <= (A - 1) * n) {
            v = u;
            break;
        }
    }
    if (v == NetworKit::none) {
        high_low_case(T, G, real_index);
        return T.getRoot();
    }
    auto V = T.add(NodeType::LEAF, real_index[v]);
    pointer[real_index[v]] = V;
    T.setRoot(V);
    std::vector<bool> is_neighbour(n);
    std::vector<int> not_neighbours;
    for (auto u : G.neighborRange(v)) {
        is_neighbour[u] = true;
    }
    for (std::size_t i = 0; i < n; i++) {
        if (i == v) {
            continue;
        }
        if (!is_neighbour[i]) {
            not_neighbours.push_back(i);
        }
    }
    std::vector<int> component(n, -1);
    std::vector<bool> is_in_vec(n, true);
    for (std::size_t i = 0; i < n; i++) {
        if (is_neighbour[i] || i == v) {
            is_in_vec[i] = false;
        }
    }
    auto components = compute_connected_components(not_neighbours, component, is_in_vec, G);
    auto gamma = compute_gamma(is_in_vec, G, component);
    components = compute_components_sorted(n, components, gamma);
    recompute_component(components, component);
    std::vector<bool> is_in_new_vec = is_in_vec;
    is_in_new_vec[v] = true;
    auto gamma_difference =
            compute_gamma_difference(components, component, gamma, is_in_vec, is_in_new_vec);
    std::vector<int> fake_index(n, -1);
    for (std::size_t i = 0; i <= components.size(); i++) {
        add(1, T, gamma_difference[i], fake_index, G, real_index);
        if (i == components.size()) {
            break;
        }
        add(0, T, components[i], fake_index, G, real_index);
    }
    return T.getRoot();
}

int current_time = 0;
std::vector<NetworKit::count> dfs_list;
NetworKit::count number_of_edges_according_to_cotree = 0;
int maximum_depth = 0;

void check_cotree_recursive(Cotree &T, NetworKit::count v, int depth) {
    maximum_depth = std::max(maximum_depth, depth);
    dfs_list.push_back(v);
    T.getNode(v).time_in = current_time++;
    if (T.getNode(v).type == NodeType::LEAF) {
        T.getNode(v).number_of_vertices_in_subtree = 1;
        T.getNode(v).time_out = current_time++;
        return;
    }
    auto child = T.getNode(v).first_child;
    int sum = 0;
    T.getNode(v).number_of_vertices_in_subtree = 0;
    while (child != NetworKit::none) {
        check_cotree_recursive(T, child, depth + 1);
        if (T.getNode(v).number == 1) {
            number_of_edges_according_to_cotree +=
                T.getNode(child).number_of_vertices_in_subtree * sum;
        }
        sum += T.getNode(child).number_of_vertices_in_subtree;
        T.getNode(v).number_of_vertices_in_subtree +=
            T.getNode(child).number_of_vertices_in_subtree;
        child = T.getNode(child).next_sibling;
    }
    T.getNode(v).time_out = current_time++;
}

bool DahlhausCographRecognition::check_cotree(Cotree &T) {
    current_time = 0;
    dfs_list.clear();
    number_of_edges_according_to_cotree = 0;
    maximum_depth = 0;
    check_cotree_recursive(T, T.getRoot(), 0);
    for (auto u : dfs_list) {
        if (T.getNode(u).parent == NetworKit::none) {
            T.getNode(u).get_up[0] = u;
        } else {
            T.getNode(u).get_up[0] = T.getNode(u).parent;
        }
    }
    int i = 1;
    for (int number = 2; number <= maximum_depth; i++, number <<= 1) {
        for (auto u : dfs_list) {
            T.getNode(u).get_up[i] = T.getNode(T.getNode(u).get_up[i - 1]).get_up[i - 1];
        }
    }
    for (auto [u, v] : graph.edgeRange()) {
        auto ancestor = lca(T, pointer[u], pointer[v], i);
        if (T.getNode(ancestor).type == NodeType::LEAF || T.getNode(ancestor).number != 1) {
            return false;
        }
    }
    return number_of_edges_according_to_cotree == graph.numberOfEdges();
}

}  // namespace Koala
