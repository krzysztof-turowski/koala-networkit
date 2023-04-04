/*
 * MaximumFlow.cpp
 *
 *  Created on: 29.03.2023
 *      Author: Michał Stobierski
 *      Ported by: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include <iostream>
#include <ranges>

#include <flow/MaximumFlow.hpp>

using edge = std::pair<NetworKit::node, NetworKit::node>;

edge reverse(const edge &p) {
    return std::make_pair(p.second, p.first);
}

namespace Koala {

MaximumFlow::MaximumFlow(NetworKit::Graph &graph, NetworKit::node s, NetworKit::node t)
    : graph(std::make_optional(graph)), source(s), target(t) { }

int MaximumFlow::getFlowSize() const {
    assureFinished();
    return flow_size;
}

int KingRaoTarjanMaximumFlow::get_visible_excess(NetworKit::node v) {
    return std::max(0, excess[v] - hidden_excess[v]);
}

NetworKit::node KingRaoTarjanMaximumFlow::get_positive_excess_node() {
    positive_excess.erase(source), positive_excess.erase(target);
    if (positive_excess.empty()) {
        return NetworKit::none;
    }
    return *positive_excess.begin();
}

void KingRaoTarjanMaximumFlow::update_positive_excess(NetworKit::node v) {
    if (get_visible_excess(v) > 0) {
        positive_excess.insert(v);
    } else {
        positive_excess.erase(v);
    }
}

int KingRaoTarjanMaximumFlow::get_flow(const edge &e) {
    auto first_parent = dynamic_tree.find_parent(e.first);
    auto second_parent = dynamic_tree.find_parent(e.second);
    if (first_parent == e.second) {
        // e in Ef
        return capacity[e] - dynamic_tree.get_value(e.first);
    }
    if (second_parent == e.first) {
        // reverse(e) in Ef
        return -(capacity[reverse(e)] - dynamic_tree.get_value(e.second));
    }
    // e in E*
    return flow[e];
}

void KingRaoTarjanMaximumFlow::set_flow(const edge &e, int c) {
    flow[e] = c, flow[reverse(e)] = -c;
}

void KingRaoTarjanMaximumFlow::saturate(const edge &e) {
    set_flow(e, capacity[e]);
    excess[e.first] -= capacity[e], excess[e.second] += capacity[e];
    update_positive_excess(e.first), update_positive_excess(e.second);
    edge_designator.response_adversary(e.first, d[e.first], e.second, d[e.first] - 1);
}

void KingRaoTarjanMaximumFlow::add_edge(const edge &e) {
    const auto &e_rev = reverse(e);
    E_star.insert(e), E_star.insert(e_rev);
    hidden_excess[e.first] -= capacity[e], hidden_excess[e.second] -= capacity[e_rev];
    update_positive_excess(e.first), update_positive_excess(e.second);
    int &d_first = d[e.first], &d_second = d[e.second];
    if (d_first > d_second) {
        saturate(e);
    } else if (d_second > d_first) {
        saturate(e_rev);
    }
}

void KingRaoTarjanMaximumFlow::cut(const edge &e) {
    set_flow(e, capacity[e] - dynamic_tree.get_value(e.first));
    dynamic_tree.cut(e.first, e.second);
    edge_designator.response_adversary(e.first, d[e.first], e.second, d[e.second]);
}

void KingRaoTarjanMaximumFlow::initialize() {
    dynamic_tree.initialize(graph->numberOfNodes());
    edge_designator.initialize(graph);
    graph->forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) {
        capacity[std::make_pair(u, v)] = w;
        hidden_excess[u] += w;
        update_positive_excess(u);
    });
    graph->forNodes([&](NetworKit::node v) {
        d[v] = 0;
    });
    for (NetworKit::count i = 0; i < graph->numberOfNodes(); i++) {
        relabel(source);
    }
    graph->forNeighborsOf(source, [&](NetworKit::node v) {
        add_edge(std::make_pair(source, v));
    });
}

std::vector<edge> KingRaoTarjanMaximumFlow::get_edges_list() {
    std::map<edge, int> edges_cost;
    graph->forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) {
        if (u == source || v == source) {
            return;
        }
        edges_cost[std::make_pair(std::min(u, v), std::max(u, v))] += w;
    });
    auto edges = std::views::keys(edges_cost);
    std::vector<edge> sorted_edges{edges.begin(), edges.end()};
    sort(sorted_edges.begin(), sorted_edges.end(), [&](const edge &a, const edge &b) {
        return edges_cost[a] < edges_cost[b];
    });
    return sorted_edges;
}

void KingRaoTarjanMaximumFlow::tree_push(NetworKit::node v, NetworKit::node u) {
    // Add current edge to the tree if needed
    if (dynamic_tree.find_root(v) == v) {
        auto e = std::make_pair(v, u);
        dynamic_tree.link(v, u, capacity[e] - get_flow(e));
    }

    // Push flow
    auto delta = std::min(
        dynamic_tree.get_minimum_path_residue_capacity(v), get_visible_excess(v));
    dynamic_tree.add_value(v, -delta);

    // Update excess
    auto path_end = dynamic_tree.find_root(v);
    excess[v] -= delta, excess[path_end] += delta;
    update_positive_excess(v), update_positive_excess(path_end);

    // Remove saturated edges
    for (auto e = dynamic_tree.find_saturated_edge(v); e.first != e.second;
            e = dynamic_tree.find_saturated_edge(e.second)) {
        cut(e);
    }
}

void KingRaoTarjanMaximumFlow::relabel(NetworKit::node v) {
    for (const auto &c : dynamic_tree.find_children(v)) {
        cut(std::make_pair(c, v));
    }

    // relabel => remove node <v,d[v]> in the game
    edge_designator.response_adversary(v, d[v]);
    d[v]++;

    // remove ineligible edges incident to v in the game
    graph->forNeighborsOf(v, [&](NetworKit::node w) {
        auto e = std::make_pair(v, w);
        // eligibility constraint
        if (!(capacity[e] - get_flow(e) > 0 && d[v] == d[w] + 1 && E_star.count(e))) {
            edge_designator.response_adversary(v, d[v], w, d[v] - 1);
        }
    });
}

void KingRaoTarjanMaximumFlow::run() {
    initialize();
    auto L = get_edges_list();
    while (!L.empty()) {
        auto e = L.back();
        L.pop_back();
        add_edge(e);
        for (NetworKit::node v = get_positive_excess_node(); v != NetworKit::none;
                v = get_positive_excess_node()) {
            auto u = edge_designator.current_edge(v, d[v]);
            if (u != NetworKit::none) {
                tree_push(v, u);
            } else {
                relabel(v);
            }
        }
    }
    flow_size = get_visible_excess(target);
    hasRun = true;
}

} /* namespace Koala */
