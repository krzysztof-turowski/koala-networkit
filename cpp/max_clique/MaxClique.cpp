#include <graph/GraphTools.hpp>

#include "max_clique/MaxClique.hpp"
#include "recognition/CographAlg.hpp"

namespace Koala {
    NetworKit::count MaxClique::recurse_run(NetworKit::count n, NetworKit::count v) {
        NetworKit::count i;
        if (cotree->nodes[v].type == 2) {
            return 1;
        } else {
            NetworKit::count l = 0, r = 0;
            if (cotree->nodes[v].left_son != -1) {
                l = recurse_run(n, cotree->nodes[v].left_son);
            }

            if (cotree->nodes[v].right_son != -1) {
                r = recurse_run(n, cotree->nodes[v].right_son);
            }

            if (cotree->nodes[v].type == 0) {
                return std::max(l, r);
            } else {
                return l + r;
            }
        }
    }

    void MaxClique::run() {
        if (cotree->prepared == false) {
            cotree->BuildTree();
        }
        size = recurse_run(n, n);
    }

    NetworKit::count MaxClique::BruteForceCliqueSize(NetworKit::Graph &Graph) {
        NetworKit::count n = Graph.numberOfNodes(), mask, i, j, ans = 0, flag;
        std::vector<NetworKit::count> st(31), clique_nodes;
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
                    if (Graph.hasEdge(clique_nodes[i], clique_nodes[j]) == false) {
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