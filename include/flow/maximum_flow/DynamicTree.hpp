/*
 * DynamicTree.hpp
 *
 *  Created on: 30.03.2023
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

#include <unordered_set>
#include <vector>

#include <networkit/graph/Graph.hpp>

#include "dynamic_tree/dyn_tree.h"

namespace Koala {

class DynamicTree {
public:
    void initialize(NetworKit::count);

    NetworKit::node find_parent(NetworKit::node);
    NetworKit::node find_root(NetworKit::node);
    std::unordered_set<NetworKit::node> find_children(NetworKit::node);

    int get_value(NetworKit::node);
    void add_value(NetworKit::node, int);

    void link(NetworKit::node, NetworKit::node, int);
    void cut(NetworKit::node, NetworKit::node);

    std::pair<NetworKit::node, NetworKit::node> find_saturated_edge(NetworKit::node);
    int get_minimum_path_residue_capacity(NetworKit::node);

private:
    std::vector<dyn_item> nodes;
    std::unordered_map<dyn_item*, NetworKit::node> ids;
    std::vector<std::unordered_set<NetworKit::node>> children;
};

} /* namespace Koala */
