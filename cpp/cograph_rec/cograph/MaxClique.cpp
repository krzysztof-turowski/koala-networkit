#include <graph/GraphTools.hpp>

#include "cograph_rec/Cotree.hpp"

namespace Koala {
    long long Cotree::MaxCliqueSize(long long n, long long v) {
        long long i;
        if (type[v] == 2) {
            return 1;
        } else {
            long long l = 0, r = 0;
            if (left_son[v] != -1) {
                l = MaxCliqueSize(n, left_son[v]);
            }

            if (right_son[v] != -1) {
                r = MaxCliqueSize(n, right_son[v]);
            }

            if (type[v] == 0) {
                return std::max(l, r);
            } else {
                return l + r;
            }
        }
    }

    long long Cotree::GetMaxCliqueSize() {
        long long n = graph->numberOfNodes();
        if (prepared == false) {
            BuildTree();
        }
        return MaxCliqueSize(n, n);
    }

    long long Cotree::BruteForceCliqueSize() {
        long long n = graph->numberOfNodes(), mask, i, j, ans = 0, flag;
        std::vector<long long> st(31), clique_nodes;
        st[0] = 1;
        for (i = 1; i <= 30; i++) {
            st[i] = st[i - 1] * 2;
        }
        for (mask = 1; mask < st[n]; mask++) {
            clique_nodes.clear();
            for (i = 0; i <= 30; i++) {
                if ((st[i] & mask) > 0) {
                    clique_nodes.push_back(i);
                }
            }
            flag = 0;
            for (i = 0; i < clique_nodes.size(); i++) {
                for (j = i + 1; j < clique_nodes.size(); j++) {
                    if (graph->hasEdge(clique_nodes[i], clique_nodes[j]) == false) {
                        flag = 1;
                        break;
                    }
                }
                if (flag == 1) {
                    break;
                }
            }
            if (flag == 0 && clique_nodes.size() > ans) {
                ans = clique_nodes.size();
            }
        }
        return ans;
    }
}