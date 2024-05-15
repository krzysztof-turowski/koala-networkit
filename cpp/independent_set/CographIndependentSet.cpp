#include <graph/GraphTools.hpp>

#include "independent_set/CographIndependentSet.hpp"

namespace Koala {
    NetworKit::count CographIndependentSet::recurse_run(NetworKit::count n, NetworKit::count v) {
        NetworKit::count i;
        if (cotree->nodes[v].type == 2) {
            return 1;
        } else {
            NetworKit::count l = 0, r = 0;
            if (cotree->nodes[v].left_son != -1) {
                l = recurse_run(n,cotree->nodes[v].left_son);
            }

            if (cotree->nodes[v].right_son != -1) {
                r = recurse_run(n,cotree->nodes[v].right_son);
            }

            if (cotree->nodes[v].type == 1) {
                return std::max(l, r);
            } else {
                return l + r;
            }
        }
    }

    void CographIndependentSet::run() {

        if (cotree->prepared == false) {
            cotree->BuildTree();
        }
        independet_set_size = recurse_run(n, n);
    }

    NetworKit::count CographIndependentSet::BruteForceIndependetSetSize(NetworKit::Graph &Graph) {
        NetworKit::count n = Graph.numberOfNodes(), mask, i, j, ans = 0, flag;
        std::vector<NetworKit::count> st(31), independet_set_nodes;
        st[0] = 1;
        for (i = 1; i <= 30; i++) {
            st[i] = st[i - 1] * 2;
        }
        for (mask = 1; mask < st[n]; mask++) {
            independet_set_nodes.clear();
            for (i = 0; i <= 30; i++) {
                if ((st[i] & mask) > 0) {
                    independet_set_nodes.push_back(i);
                }
            }
            flag = 0;
            for (i = 0; i < independet_set_nodes.size(); i++) {
                for (j = i + 1; j < independet_set_nodes.size(); j++) {
                    if (Graph.hasEdge(independet_set_nodes[i], independet_set_nodes[j])) {
                        flag = 1;
                        break;
                    }
                }
                if (flag == 1) {
                    break;
                }
            }
            if (flag == 0 && independet_set_nodes.size() > ans) {
                ans = independet_set_nodes.size();
            }
        }
        return ans;
    }
}