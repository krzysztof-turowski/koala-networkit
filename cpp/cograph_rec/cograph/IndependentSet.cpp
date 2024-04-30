#include <graph/GraphTools.hpp>

#include "cograph_rec/Cotree.hpp"

namespace Koala {
    long long Cotree::MaxIndependetSetSize(long long n, long long v) {
        long long i;
        if (type[v] == 2) {
            return 1;
        } else {
            long long l = 0, r = 0;

            if (left_son[v] != -1) {
                l = MaxIndependetSetSize(n, left_son[v]);
            }

            if (right_son[v] != -1) {
                r = MaxIndependetSetSize(n, right_son[v]);
            }

            if (type[v] == 1) {
                return std::max(l, r);
            } else {
                return l + r;
            }
        }
    }

    long long Cotree::GetMaxIndependetSetSize() {
        long long n = graph->numberOfNodes();
        if (prepared == false) {
            BuildTree();
        }
        return MaxIndependetSetSize(n, n);
    }

    long long Cotree::BruteForceIndependetSetSize() {
        long long n = graph->numberOfNodes(), mask, i, j, ans = 0, flag;
        std::vector<long long> st(31), independet_set_nodes;
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
                    if (graph->hasEdge(independet_set_nodes[i], independet_set_nodes[j])) {
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