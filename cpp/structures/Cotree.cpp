#include <graph/GraphTools.hpp>

#include "structures/Cotree.hpp"

namespace Koala {
    Cotree::Cotree(NetworKit::Graph &Graph) {
        graph = &Graph;
        prepared = false;
    }

    void Cotree::BuildTree() {
        reverse(order.begin(), order.end());
        NetworKit::count n = graph->numberOfNodes(), i;
        nodes.resize(2 * n, CoNode(0, 0, 0));

        nodes[n].left_son = order[0].first.first;
        nodes[n].right_son = -1;
        nodes[n].parent = -1;
        nodes[n].type = 0;

        nodes[order[0].first.first].left_son = -1;
        nodes[order[0].first.first].right_son = -1;
        nodes[order[0].first.first].type = 2;
        nodes[order[0].first.first].parent = n;

        for (i = 1; i < n; i++) {
            nodes[order[i].first.first].left_son = -1;
            nodes[order[i].first.first].right_son = -1;
            nodes[order[i].first.first].type = 2;
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
            nodes[n + i].type = order[i].second;
        }
        prepared = true;
    }
}