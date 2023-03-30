/*
 * DynamicTree.hpp
 *
 *  Created on: 30.03.2023
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

#include <map>
#include <vector>

#include "dynamic_trees/dyn_tree.h"

namespace Koala {

class DynamicTree {
public:
    void initialize(NetworKit::count n) {
        nodes.resize(n + 1);
        children.resize(n + 1);
        for (int v = 0; v <= n; ++v) {
            ids[&nodes[v]] = v;
            dyn_make_tree(&nodes[v], 0);
        }
    }

    NetworKit::node find_parent(NetworKit::node v) {
        return ids[dyn_find_father(&nodes[v + 1])] - 1;
    }

    NetworKit::node find_root(NetworKit::node v) {
        return ids[dyn_find_root(&nodes[v + 1])] - 1;
    }

    std::unordered_set<NetworKit::node> find_children(NetworKit::node v) {
        return children[v + 1];
    }

    int get_value(NetworKit::node v) {
        return dyn_find_value(&nodes[v + 1]);
    }

    void add_value(NetworKit::node v, int delta) {
        dyn_add_value(&nodes[v + 1], delta);
    }

    void link(NetworKit::node u, NetworKit::node v, int capacity) {
        dyn_link(&nodes[u + 1], &nodes[v + 1], capacity);
        children[v + 1].insert(u);
    }

    void cut(NetworKit::node u, NetworKit::node v) {
        dyn_cut(&nodes[u + 1]);
        children[v + 1].erase(u);
    }

    std::pair<NetworKit::node, NetworKit::node> find_saturated_edge(NetworKit::node v) {
        dyn_item *start = &nodes[v + 1];
        dyn_item *bottleneck = dyn_find_bottleneck(start, 0), *root = dyn_find_root(start);
        if (bottleneck != root) {
            return std::make_pair(ids[bottleneck] - 1, ids[dyn_find_father(bottleneck)] - 1);
        }
        return std::make_pair(NetworKit::none, NetworKit::none);
    }

    int get_minimum_path_residue_capacity(NetworKit::node v) {
        v += 1;
        int p = 0, q = std::numeric_limits<int>::max() - 1;
        dyn_item *root = dyn_find_root(&nodes[v]);
        while (p < q) {
            int mid = (p + q) / 2;
            dyn_item *bot = dyn_find_bottleneck(&nodes[v], mid);
            if (bot == root) {
                p = mid + 1;
            } else {
                q = mid;
            }
        }
        return p;
    }

private:
    std::vector<dyn_item> nodes;
    std::unordered_map<dyn_item*, NetworKit::node> ids;
    std::vector<std::unordered_set<NetworKit::node>> children;
};

} /* namespace Koala */
