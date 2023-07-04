/*
 * DynamicTree.cpp
 *
 *  Created on: 30.03.2023
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#include <flow/maximum_flow/DynamicTree.hpp>

#include "dynamic_tree/dyn_tree.h"

namespace Koala {

void DynamicTree::initialize(NetworKit::count n) {
    nodes.resize(n), children.resize(n);
    for (NetworKit::count v = 0; v < n; ++v) {
        ids[&nodes[v]] = v;
        dyn_make_tree(&nodes[v], 0);
    }
    ids[NULL] = NetworKit::none;
}

NetworKit::node DynamicTree::find_parent(NetworKit::node v) {
    return ids[dyn_find_father(&nodes[v])];
}

NetworKit::node DynamicTree::find_root(NetworKit::node v) {
    return ids[dyn_find_root(&nodes[v])];
}

std::unordered_set<NetworKit::node> DynamicTree::find_children(NetworKit::node v) {
    return children[v];
}

int DynamicTree::get_value(NetworKit::node v) {
    return dyn_find_value(&nodes[v]);
}

void DynamicTree::add_value(NetworKit::node v, int delta) {
    dyn_add_value(&nodes[v], delta);
}

void DynamicTree::link(NetworKit::node u, NetworKit::node v, int capacity) {
    dyn_link(&nodes[u], &nodes[v], capacity);
    children[v].insert(u);
}

void DynamicTree::cut(NetworKit::node u, NetworKit::node v) {
    dyn_cut(&nodes[u]);
    children[v].erase(u);
}

std::pair<NetworKit::node, NetworKit::node> DynamicTree::find_saturated_edge(NetworKit::node v) {
    dyn_item *start = &nodes[v];
    dyn_item *bottleneck = dyn_find_bottleneck(start, 0), *root = dyn_find_root(start);
    if (bottleneck != root) {
        return std::make_pair(ids[bottleneck], ids[dyn_find_father(bottleneck)]);
    }
    return std::make_pair(NetworKit::none, NetworKit::none);
}

int DynamicTree::get_minimum_path_residue_capacity(NetworKit::node v) {
    int p = 0, q = std::numeric_limits<int>::max() - 1;
    dyn_item *root = dyn_find_root(&nodes[v]);
    while (p < q) {
        int mid = (p + q) / 2;
        dyn_item *bottleneck = dyn_find_bottleneck(&nodes[v], mid);
        if (bottleneck == root) {
            p = mid + 1;
        } else {
            q = mid;
        }
    }
    return p;
}

} /* namespace Koala */
