#include <graph/GraphTools.hpp>

#include "independent_set/CographIndependentSet.hpp"

namespace Koala {
    NetworKit::count CographIndependentSet::recurse_run(NetworKit::count v)  {
        CoNode &V = cotree.getNode(v);
        if (V.type == NodeType::LEAF) {
            independent_set_size[v] = 1;
        } else {
            NetworKit::count  l=0, r=0;
            if (V.left_son != NetworKit::none) {
                l = recurse_run(V.left_son);
            }

            if (V.right_son != NetworKit::none) {
                r = recurse_run(V.right_son);
            }

            if (V.type == NodeType::COMPLEMENT_NODE) {
                independent_set_size[v] = std::max(l,r);
            } else {
                independent_set_size[v] = l + r;
            }
        }
        return independent_set_size[v];
    }

    void CographIndependentSet::add_to_set(NetworKit::count v) {
        if(v == NetworKit::none)
        {
            return;
        }
        CoNode &V = cotree.getNode(v);
        if (V.type == NodeType::LEAF) {
           independentSet.insert(v);
        } else {
            NetworKit::count  l=0, r=0;
            if (V.left_son != NetworKit::none) {
                l = independent_set_size[V.left_son];
            }

            if (V.right_son != NetworKit::none) {
                r = independent_set_size[V.right_son];
            }

            if (V.type == NodeType::COMPLEMENT_NODE) {
                if (l >= r) {
                    add_to_set(V.left_son);
                } else {
                    add_to_set(V.right_son);
                }
            } else {
                add_to_set(V.left_son);
                add_to_set(V.right_son);
            }
        }
    }

    void CographIndependentSet::run() {
        hasRun = true;
        NetworKit::count n = cotree.graph->numberOfNodes();
        independent_set_size.resize(2*n+1);
        recurse_run(n);
        add_to_set(n);
    }

    NetworKit::count CographIndependentSet::bruteForceIndependetSetSize(NetworKit::Graph &Graph) {
        NetworKit::count n = Graph.numberOfNodes(), ans = 0, flag;
        std::vector<NetworKit::count> st(31), independet_set_nodes;
        st[0] = 1;
        for (int i = 1; i <= 30; i++) {
            st[i] = st[i - 1] * 2;
        }
        for (int mask = 1; mask < st[n]; mask++) {
            independet_set_nodes.clear();
            for (int i = 0; i <= 30; i++) {
                if ((st[i] & mask) > 0) {
                    independet_set_nodes.push_back(i);
                }
            }
            flag = 0;
            for (int i = 0; i < independet_set_nodes.size(); i++) {
                for (int j = i + 1; j < independet_set_nodes.size(); j++) {
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
} /* namespace Koala */
