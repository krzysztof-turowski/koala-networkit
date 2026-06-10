/*
 * DahlhausCographRecognition.cpp
 *
 *  Created on: 07.06.2024
 *      Author: fixikmila
 */

#include <algorithm>
#include <optional>
#include <vector>

#include <graph/GraphTools.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <recognition/CographRecognition.hpp>
#include <structures/Cotree.hpp>
#include <structures/LCA.hpp>

namespace Koala {

void DahlhausCographRecognition::run() {
    is_cograph = State::COGRAPH;
    hasRun = true;
    T.clear();
    T.reserve(graph.numberOfNodes() * 4);
    covertex.assign(graph.upperNodeIdBound(), NetworKit::none);
    auto root = build_cotree(graph);
    T.setRoot(root);
    if (is_cograph == State::COGRAPH) {
        is_cograph = check_cotree() ? State::COGRAPH : State::NOT_COGRAPH;
    }
}

inline void DahlhausCographRecognition::add(
        NodeType node_type, std::vector<NetworKit::node> &vec, NetworKit::Graph &G) {
    auto u2 = T.add(node_type);
    T.addChild(u2, T.getRoot());
    T.setRoot(u2);
    if (vec.empty()) {
        return;
    }
    auto C = NetworKit::GraphTools::subgraphFromNodes(G, vec.begin(), vec.end());
    auto root = T.getRoot();
    auto subtree_root = build_cotree(C);
    if (is_cograph != State::COGRAPH) {
        return;
    }
    T.setRoot(root);
    if (node_type == NodeType::UNION_NODE) {
        T.addChild(u2, subtree_root);
    } else {
        auto u1 = T.add(NodeType::UNION_NODE);
        T.addChild(u2, u1);
        T.addChild(u1, subtree_root);
    }
    T.setRoot(u2);
}

void recompute_component(
        std::vector<std::vector<NetworKit::node>> &components,
        std::vector<NetworKit::count> &component) {
    for (std::size_t i = 0; i < components.size(); i++) {
        for (std::size_t j = 0; j < components[i].size(); j++) {
            component[components[i][j]] = i;
        }
    }
}

std::vector<std::vector<NetworKit::node>> compute_gamma(
        std::vector<bool> &is_in_vec, std::vector<std::vector<NetworKit::node>> &components,
        NetworKit::Graph &G) {
    std::vector<NetworKit::count> component(G.upperNodeIdBound(), NetworKit::none);
    recompute_component(components, component);
    std::vector<std::vector<NetworKit::node>> gamma(components.size());
    std::vector<NetworKit::count> last_component(G.upperNodeIdBound(), NetworKit::none);
    for (std::size_t i = 0; i < components.size(); i++) {
        for (auto v : components[i]) {
            for (auto u : G.neighborRange(v)) {
                if (is_in_vec[u] && component[u] == i) {
                    continue;
                }
                if (last_component[u] == i) {
                    continue;
                }
                last_component[u] = i;
                gamma[i].push_back(u);
            }
        }
    }
    return gamma;
}

std::vector<std::vector<NetworKit::node>> compute_components_sorted(
        NetworKit::count n, std::vector<std::vector<NetworKit::node>> &components,
        std::vector<std::vector<NetworKit::node>> &gamma) {
    std::vector<std::vector<NetworKit::count>> count_sort(n);
    std::vector<std::vector<NetworKit::node>> components_sorted;
    for (std::size_t i = 0; i < components.size(); i++) {
        count_sort[gamma[i].size()].push_back(i);
    }
    for (auto it = count_sort.rbegin(); it != count_sort.rend(); ++it) {
        for (auto value : *it) {
            components_sorted.push_back(components[value]);
        }
    }
    return components_sorted;
}

std::optional<std::vector<std::vector<NetworKit::node>>>
compute_gamma_difference(
        std::vector<std::vector<NetworKit::node>> &components,
        std::vector<std::vector<NetworKit::node>> &gamma, NetworKit::Graph &G,
        std::vector<bool> &is_in_new_vec) {
    std::vector<std::vector<NetworKit::node>> gamma_difference(components.size() + 1);
    std::vector<NetworKit::node> last_position_where_met(
        G.upperNodeIdBound(), NetworKit::none);
    for (std::size_t i = 0; i < components.size(); i++) {
        for (auto a : gamma[i]) {
            if (last_position_where_met[a] != NetworKit::none
                    && last_position_where_met[a] + 1 != i) {
                return std::nullopt;
            }
            last_position_where_met[a] = i;
        }
    }
    for (auto i : G.nodeRange()) {
        if (is_in_new_vec[i]) {
            continue;
        }
        const auto position = last_position_where_met[i] == NetworKit::none
                ? 0 : last_position_where_met[i] + 1;
        gamma_difference[position].push_back(i);
    }
    return gamma_difference;
}

void reverse_cotree(Cotree &T, NetworKit::node v) {
    if (T.getNode(v).type != NodeType::LEAF) {
        T.getNode(v).type = T.getNode(v).type == NodeType::UNION_NODE
            ? NodeType::COMPLEMENT_NODE : NodeType::UNION_NODE;
    }
    auto u = T.getNode(v).first_child;
    while (u != NetworKit::none) {
        reverse_cotree(T, u);
        u = T.getNode(u).next_sibling;
    }
}

void DahlhausCographRecognition::big_component(
        NetworKit::Graph &G, std::vector<NetworKit::node> &vec) {
    NetworKit::Graph GC = Koala::GraphTools::toComplement(G);
    NetworKit::count n = GC.numberOfNodes();
    std::vector<std::vector<NetworKit::node>> components;
    if (!vec.empty()) {
        auto induced = NetworKit::GraphTools::subgraphFromNodes(GC, vec.begin(), vec.end());
        NetworKit::ConnectedComponents connected_components(induced);
        connected_components.run();
        components = connected_components.getComponents();
    }
    for (auto c : components) {
        if (c.size() * A > 2 * n + A) {
            is_cograph = State::NOT_COGRAPH;
            return;
        }
        auto GI = NetworKit::GraphTools::subgraphFromNodes(GC, c.begin(), c.end());
        auto root = T.getRoot();
        auto subtree_root = build_cotree(GI);
        if (is_cograph != State::COGRAPH) {
            return;
        }
        T.setRoot(root);
        reverse_cotree(T, subtree_root);  // reverse TI
        T.addChild(T.getRoot(), subtree_root);
    }
}

void DahlhausCographRecognition::high_low_case(NetworKit::Graph &G) {
    if (is_cograph != State::COGRAPH) {
        return;
    }
    NetworKit::count n = G.numberOfNodes();
    std::vector<NetworKit::count> degree(G.upperNodeIdBound());
    for (auto u : G.nodeRange()) {
        degree[u] = G.degree(u);
    }
    auto V = T.add(NodeType::UNION_NODE);
    T.setRoot(V);
    // compute low components
    // sort and compute gamma difference
    // 0-low components 1-high components
    // if high component is big, then call big_components
    std::vector<bool> low(G.upperNodeIdBound());
    std::vector<NetworKit::node> vec;
    for (auto u : G.nodeRange()) {
        if (degree[u] * A <= n) {
            low[u] = true;
            vec.push_back(u);
        }
    }
    std::vector<std::vector<NetworKit::node>> components;
    if (!vec.empty()) {
        auto induced = NetworKit::GraphTools::subgraphFromNodes(G, vec.begin(), vec.end());
        NetworKit::ConnectedComponents connected_components(induced);
        connected_components.run();
        components = connected_components.getComponents();
    }
    auto gamma = compute_gamma(low, components, G);
    components = compute_components_sorted(n, components, gamma);
    gamma = compute_gamma(low, components, G);
    auto gamma_difference = compute_gamma_difference(components, gamma, G, low);
    if (!gamma_difference) {
        is_cograph = State::NOT_COGRAPH;
        return;
    }
    auto &gamma_difference_value = *gamma_difference;

    for (std::size_t i = 0; i <= components.size(); i++) {
        bool special_case_big_component = true;
        if (gamma_difference_value[i].size() * A > (A - 1) * n) {
            std::vector<bool> is_in_gamma_difference(G.upperNodeIdBound());
            for (auto u : gamma_difference_value[i]) {
                is_in_gamma_difference[u] = true;
            }
            for (auto u : gamma_difference_value[i]) {
                NetworKit::count sum = 0;
                for (auto v : G.neighborRange(u)) {
                    if (!is_in_gamma_difference[v]) {
                        continue;
                    }
                    sum++;
                }
                sum = gamma_difference_value[i].size() - 1 - sum;
                if (sum * A >= n) {
                    special_case_big_component = false;
                    break;
                }
            }
        } else {
            special_case_big_component = false;
        }
        if (special_case_big_component) {
            std::vector<NetworKit::node> empty;
            add(NodeType::COMPLEMENT_NODE, empty, G);
            big_component(G, gamma_difference_value[i]);
            if (is_cograph != State::COGRAPH) {
                return;
            }
        } else {
            add(NodeType::COMPLEMENT_NODE, gamma_difference_value[i], G);
        }
        if (i == components.size()) {
            break;
        }
        if (components[i].size() * A > 2 * n + A) {
            is_cograph = State::NOT_COGRAPH;
            return;
        }
        add(NodeType::UNION_NODE, components[i], G);
    }
}

NetworKit::node DahlhausCographRecognition::build_cotree(NetworKit::Graph &G) {
    NetworKit::count n = G.numberOfNodes();
    T.reserve(3 * n);
    if (n == 1) {
        auto v = *G.nodeRange().begin();
        auto cv = covertex[v] = T.add(NodeType::LEAF, v);
        T.setRoot(cv);
        return cv;
    }
    if (is_cograph != State::COGRAPH) {
        return T.getRoot();
    }

    auto node_range = G.nodeRange();
    auto found = std::find_if(node_range.begin(), node_range.end(), [&](NetworKit::node u) {
        auto size = G.degree(u);
        return A * size >= n && size * A <= (A - 1) * n;
    });
    auto v = found != node_range.end() ? *found : NetworKit::none;
    if (v == NetworKit::none) {
        high_low_case(G);
        return T.getRoot();
    }
    auto V = T.add(NodeType::LEAF, v);
    covertex[v] = V;
    T.setRoot(V);
    std::vector<bool> is_neighbour(G.upperNodeIdBound());
    std::vector<NetworKit::node> not_neighbours;
    for (auto u : G.neighborRange(v)) {
        is_neighbour[u] = true;
    }
    std::vector<bool> is_in_vec(G.upperNodeIdBound());
    for (auto i : G.nodeRange()) {
        if (i != v && !is_neighbour[i]) {
            not_neighbours.push_back(i);
            is_in_vec[i] = true;
        }
    }
    std::vector<std::vector<NetworKit::node>> components;
    if (!not_neighbours.empty()) {
        auto induced = NetworKit::GraphTools::subgraphFromNodes(
            G, not_neighbours.begin(), not_neighbours.end());
        NetworKit::ConnectedComponents connected_components(induced);
        connected_components.run();
        components = connected_components.getComponents();
    }
    auto gamma = compute_gamma(is_in_vec, components, G);
    components = compute_components_sorted(n, components, gamma);
    gamma = compute_gamma(is_in_vec, components, G);
    is_in_vec[v] = true;
    auto gamma_difference = compute_gamma_difference(components, gamma, G, is_in_vec);
    if (!gamma_difference) {
        is_cograph = State::NOT_COGRAPH;
        return T.getRoot();
    }
    for (std::size_t i = 0; i <= components.size(); i++) {
        add(NodeType::COMPLEMENT_NODE, (*gamma_difference)[i], G);
        if (i == components.size()) {
            break;
        }
        add(NodeType::UNION_NODE, components[i], G);
    }
    return T.getRoot();
}

NetworKit::count count_edges(Cotree &T) {
    NetworKit::count total_edges = 0;
    std::vector<NetworKit::count> vertices(T.upperNodeIdBound(), 0);
    std::vector<NetworKit::node> order, stack = {T.getRoot()};
    while (!stack.empty()) {
        auto v = stack.back();
        stack.pop_back();
        order.push_back(v);
        auto child = T.getNode(v).first_child;
        while (child != NetworKit::none) {
            stack.push_back(child);
            child = T.getNode(child).next_sibling;
        }
    }
    for (auto it = order.rbegin(); it != order.rend(); ++it) {
        auto v = *it;
        if (T.getNode(v).type == NodeType::LEAF) {
            vertices[v] = 1;
            continue;
        }
        auto child = T.getNode(v).first_child;
        NetworKit::count sum = 0;
        while (child != NetworKit::none) {
            if (T.getNode(v).type == NodeType::COMPLEMENT_NODE) {
                total_edges += vertices[child] * sum;
            }
            sum += vertices[child];
            vertices[v] += vertices[child];
            child = T.getNode(child).next_sibling;
        }
    }
    return total_edges;
}

std::vector<NetworKit::node> get_parents(Cotree &T) {
    std::vector<NetworKit::node> parent(T.upperNodeIdBound(), NetworKit::none);
    std::vector<NetworKit::node> stack = {T.getRoot()};
    while (!stack.empty()) {
        auto u = stack.back();
        stack.pop_back();
        auto child = T.getNode(u).first_child;
        while (child != NetworKit::none) {
            parent[child] = u;
            stack.push_back(child);
            child = T.getNode(child).next_sibling;
        }
    }
    return parent;
}

bool DahlhausCographRecognition::check_cotree() {
    if (count_edges(T) != graph.numberOfEdges()) {
        return false;
    }

    LCA lca(get_parents(T), T.getRoot());
    for (auto [u, v] : graph.edgeRange()) {
        if (covertex[u] == NetworKit::none || covertex[v] == NetworKit::none) {
            return false;
        }
        auto ancestor = lca.query(covertex[u], covertex[v]);
        if (T.getNode(ancestor).type != NodeType::COMPLEMENT_NODE) {
            return false;
        }
    }

    return true;
}

}  // namespace Koala
