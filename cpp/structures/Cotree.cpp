#include <graph/GraphTools.hpp>

#include "structures/Cotree.hpp"

namespace Koala {

Cotree::Cotree(NetworKit::Graph &Graph) {
    graph = &Graph;
    prepared = false;
}

void Cotree::buildTree() {
    reverse(order.begin(), order.end());
    NetworKit::count n = graph->numberOfNodes(), i;
    nodes.resize(2 * n, Conode(0, 0, 0));

    nodes[n].left_son = order[0].first.first;
    nodes[n].right_son = NetworKit::none;
    nodes[n].parent = NetworKit::none;
    nodes[n].type = NodeType::UNION_NODE;

    nodes[order[0].first.first].left_son = NetworKit::none;
    nodes[order[0].first.first].right_son = NetworKit::none;
    nodes[order[0].first.first].type = NodeType::LEAF;
    nodes[order[0].first.first].parent = n;

    for (i = 1; i < n; i++) {
        nodes[order[i].first.first].left_son = NetworKit::none;
        nodes[order[i].first.first].right_son = NetworKit::none;
        nodes[order[i].first.first].type = NodeType::LEAF;
        nodes[order[i].first.first].parent = n + i;

        NetworKit::count p = nodes[order[i].first.second].parent;
        if (nodes[p].left_son == order[i].first.second) {
            nodes[p].left_son = n + i;
        } else {
            nodes[p].right_son = n + i;
        }

        nodes[order[i].first.second].parent = n + i;

        nodes[n + i].left_son = order[i].first.first;
        nodes[n + i].right_son = order[i].first.second;
        nodes[n + i].parent = p;
        if (order[i].second == 1) {
            nodes[n + i].type = NodeType::COMPLEMENT_NODE;
        } else {
            nodes[n + i].type = NodeType::UNION_NODE;
        }
    }
    prepared = true;
}

} /* namespace Koala */
