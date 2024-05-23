#include <graph/GraphTools.hpp>

#include "max_clique/MaxClique.hpp"
#include "recognition/CographAlg.hpp"

namespace Koala {
    NetworKit::count MaxClique::recurse_run(NetworKit::count v) {
        CoNode V(0,0,0);
        V=cotree.getNode(v);
        if (V.type == NodeType::LEAF) {
            return 1;
        } else {
            NetworKit::count l = 0, r = 0;
            if (V.left_son != NetworKit::none) {
                l = recurse_run(V.left_son);
            }

            if (V.right_son != NetworKit::none) {
                r = recurse_run(V.right_son);
            }

            if (V.type == NodeType::UNION_NODE) {
                return std::max(l, r);
            } else {
                return l + r;
            }
        }
    }

    void MaxClique::run() {
        if (cotree.prepared == false) {
            cotree.buildTree();
        }
        size = recurse_run(cotree.graph->numberOfNodes());
    }

    NetworKit::count MaxClique::bruteForceCliqueSize(NetworKit::Graph &Graph) {
        NetworKit::count n = Graph.numberOfNodes(), mask, ans = 0, flag;
        std::vector<NetworKit::count> st(31), clique_nodes;
        st[0] = 1;
        for (int i = 1; i <= 30; i++) {
            st[i] = st[i - 1] * 2;
        }
        for (mask = 1; mask < st[n]; mask++) {
            clique_nodes.clear();
            for (int i = 0; i <= 30; i++) {
                if ((st[i] & mask) > 0) {
                    clique_nodes.push_back(i);
                }
            }
            flag = 0;

            for (int i = 0; i < clique_nodes.size(); i++) {
                for (int j = i + 1; j < clique_nodes.size(); j++) {
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