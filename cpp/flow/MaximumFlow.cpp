/*
 * MaximumFlow.cpp
 *
 *  Created on: 29.03.2023
 *      Author: Michał Stobierski
 *      Ported by: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include <iostream>

#include <flow/MaximumFlow.hpp>

#include <flow/maximum_flow/KrtEdgeDesignator.cpp>
#include <flow/maximum_flow/dynamic_trees/dyn_tree.cpp>
#include <flow/maximum_flow/dynamic_trees/dyn_splay.cpp>

namespace Koala {

using edge = std::pair<NetworKit::node, NetworKit::node>;

MaximumFlow::MaximumFlow(const NetworKit::Graph &graph, NetworKit::node s, NetworKit::node t)
    : graph(std::make_optional(graph)), source(s), target(t) { }

int MaximumFlow::getFlowSize() const {
    assureFinished();
    return flow_size;
}

edge MaximumFlow::reverse(const edge &p) {
    return std::make_pair(p.second, p.first);
}

int KingRaoTarjanMaximumFlow::get_visible_excess(NetworKit::node v) {
    return std::max(0, excess[v] - hidden_excess[v]);
}

void KingRaoTarjanMaximumFlow::update_positive_excess_map(NetworKit::node v) {
    if (get_visible_excess(v) > 0) {
        positive_excess.insert(v);
    } else {
        positive_excess.erase(v);
    }
}

NetworKit::node KingRaoTarjanMaximumFlow::get_positive_excess_node() {
    positive_excess.erase(source);
    positive_excess.erase(target);
    if (!positive_excess.empty()) {
        return *positive_excess.begin();
    }
    return NetworKit::none;
}

int KingRaoTarjanMaximumFlow::get_flow(const edge &e) {
    auto first_parent = dynamic_tree.find_parent(e.first);
    auto second_parent = dynamic_tree.find_parent(e.second);

    int res;
    if (first_parent == e.second) {
        // e in Ef
        res = capacity[e] - dynamic_tree.get_value(e.first);
    } else if (second_parent == e.first) {
        // reverse(e) in Ef
        res = -(capacity[reverse(e)] - dynamic_tree.get_value(e.second));
    } else {
        // e in E*
        res = flow[e];
    }
    logger.log("f", e, res);
    return res;
}

void KingRaoTarjanMaximumFlow::set_flow(const edge &e, int c) {
    logger.log("setFlow", e, c);
    flow[e] = c, flow[reverse(e)] = -c;
}

void KingRaoTarjanMaximumFlow::saturate(const edge &e) {
    logger.log("saturate", e, capacity[e]);
    set_flow(e, capacity[e]);
    excess[e.first] -= capacity[e], excess[e.second] += capacity[e];
    update_positive_excess_map(e.first), update_positive_excess_map(e.second);

    // sat push => adversary edge kill
    edge_designator.response_adversary(e.first, d[e.first], e.second, d[e.first] - 1, false);
}

void KingRaoTarjanMaximumFlow::add_edge(const edge &e) {
    const auto &e_rev = reverse(e);
    E_star.insert(e), E_star.insert(e_rev);
    hidden_excess[e.first] -= capacity[e], hidden_excess[e.second] -= capacity[e_rev];
    update_positive_excess_map(e.first), update_positive_excess_map(e.second);
    int &d_first = d[e.first], &d_second = d[e.second];
    if (d_first > d_second) {
        saturate(e);
    } else if (d_second > d_first) {
        saturate(e_rev);
    }
}

void KingRaoTarjanMaximumFlow::initialize() {
    dynamic_tree.initialize(graph->numberOfNodes());
    edge_designator.initialize(graph);
    graph->forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) {
        capacity[std::make_pair(u, v)] = w;
        hidden_excess[u] += w;
        update_positive_excess_map(u);
    });
    graph->forNodes([&](NetworKit::node v) {
        d[v] = 0;
    });
    for (int i = 0; i < graph->numberOfNodes(); i++) {
        relabel(source);
    }
    graph->forNeighborsOf(source, [&](NetworKit::node v) {
        add_edge(std::make_pair(source, v));
    });
}

std::vector<edge> KingRaoTarjanMaximumFlow::get_edges_list() {
    logger.log("initEdgesList");
    std::map<edge, int> edge_costs;
    std::set<edge> edges;
    graph->forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) {
        if (u == source || v == source) {
            return;
        }
        auto e = std::make_pair(std::min(u, v), std::max(u, v));
        edge_costs[e] += w;
        edges.insert(e);
    });
    std::vector<edge> sorted_edges(edges.begin(), edges.end());
    sort(sorted_edges.begin(), sorted_edges.end(), [&](const edge &a, const edge &b) {
        return edge_costs[a] < edge_costs[b];
    });
    return sorted_edges;
}

void KingRaoTarjanMaximumFlow::cut(const edge &e) {
    logger.log("cut", e, -1);
    set_flow(e, capacity[e] - dynamic_tree.get_value(e.first));
    dynamic_tree.cut(e.first, e.second);

    // Sat push => adversary edge kill
    edge_designator.response_adversary(e.first, d[e.first], e.second, d[e.second], false);
}

void KingRaoTarjanMaximumFlow::tree_push(NetworKit::node v, NetworKit::node u) {
    logger.log("treePush", v);
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
    update_positive_excess_map(v), update_positive_excess_map(path_end);

    // Remove saturated edges
    for (auto e = dynamic_tree.find_saturated_edge(v); e.first != e.second;
            e = dynamic_tree.find_saturated_edge(e.second)) {
        cut(e);
    }
}

void KingRaoTarjanMaximumFlow::relabel(NetworKit::node v) {
    logger.log("relabel", v, d[v] + 1);
    for (const auto &c : dynamic_tree.find_children(v)) {
        cut(std::make_pair(c, v));
    }
    d[v]++;

    // relabel => remove node <v,d[v]> in the game
    edge_designator.response_adversary(-1, -1, v, d[v] - 1, true);

    // remove ineligible edges incident to v in the game
    graph->forNeighborsOf(v, [&](NetworKit::node w) {
        auto e = std::make_pair(v, w);
        // eligibility constraint
        if (!(capacity[e] - get_flow(e) > 0 && d[v] == d[w] + 1 && E_star.count(e))) {
            edge_designator.response_adversary(v, d[v], w, d[v] - 1, false);
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
        for (int v = get_positive_excess_node(); v != NetworKit::none; v = get_positive_excess_node()) {
            auto u = edge_designator.current_edge(v, d[v]);
            logger.log("currentEdge", v, u);
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
