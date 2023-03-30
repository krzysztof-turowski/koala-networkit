/*
 * MaximumFlow.cpp
 *
 *  Created on: 29.03.2023
 *      Author: Michał Stobierski
 *      Ported by: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include <iostream>

#include <flow/MaximumFlow.hpp>

#include "krt/king_rao_tarjan.hpp"

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
    return std::max(0, excess[v] - excess_hidden[v]);
}

void KingRaoTarjanMaximumFlow::update_positive_excess_map(NetworKit::node v) {
    if (get_visible_excess(v) > 0) {
        positive_excess.insert(v);
    } else {
        positive_excess.erase(v);
    }
}

// int get_flow(const edge &e) {
    // auto first_parent = dynamic_tree.find_father(e.first);
    // auto second_parent = dynamic_tree.find_father(e.second);

    // if (first_parent == e.second) {
        // return cap[e] - dynamic_tree.find_value(e.first);
    // }
    // if (secondParent == e.first) {
        // return -(cap[rev(e)] - dynamic_tree.find_value(e.second);
    // }
    // return flow[e];
// }

void KingRaoTarjanMaximumFlow::set_flow(const edge &e, int c) {
    flow[e] = c, flow[reverse(e)] = -c;
}

void KingRaoTarjanMaximumFlow::saturate(const edge &e) {
    set_flow(e, capacity[e]);
    excess[e.first] -= capacity[e], excess[e.second] += capacity[e];
    update_positive_excess_map(e.first);
    update_positive_excess_map(e.second);

    // sat push => adversary edge kill
    // edgeDesignator->response_adversary(e.first, d[e.first], e.second, d[e.first] - 1, false);
}

void KingRaoTarjanMaximumFlow::add_edge(const edge &e) {
    const edge &e_rev = reverse(e);
    E_star.insert(e), E_star.insert(e_rev);
    excess_hidden[e.first] -= capacity[e], excess_hidden[e.second] -= capacity[e_rev];
    update_positive_excess_map(e.first);
    update_positive_excess_map(e.second);

    int &d_first = d[e.first], &d_second = d[e.second];
    if (d_first > d_second) {
        saturate(e);
    } else if (d_second > d_first) {
        saturate(e_rev);
    }
}

void KingRaoTarjanMaximumFlow::initialize() {
    graph->forEdges([&](NetworKit::node u, NetworKit::node v, NetworKit::edgeweight w) {
        std::cout << "EDGE " << u + 1 << " " << v + 1 << " with weight " << w << std::endl;
        capacity[std::make_pair(u, v)] = w;
        excess_hidden[u] += w;
        update_positive_excess_map(u);
    });
    graph->forNodes([&](NetworKit::node v) {
        d[v] = 0;
    });
    // dynamic_tree.initialize();

    for (int i = 0; i < graph->numberOfNodes(); i++) {
        // relabel(source);
    }
    graph->forNeighborsOf(source, [&](NetworKit::node v) {
        add_edge(std::make_pair(source, v));
    });
}

void KingRaoTarjanMaximumFlow::run() {
    std::cout << "KingRaoTarjanMaximumFlow::run() started" << std::endl;
    flow_size = king_rao_tarjan(graph, source, target);
    std::cout << "KingRaoTarjanMaximumFlow::run() finished" << std::endl;
    hasRun = true;
}

} /* namespace Koala */
