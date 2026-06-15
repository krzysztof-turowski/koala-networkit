#include <algorithm>
#include <vector>

#include <graph/GraphTools.hpp>
#include <independent_set/CographIndependentSet.hpp>

namespace Koala {

void CographIndependentSet::recurse_run() {
    while (!st.empty()) {
        NetworKit::node v = st.top();
        Conode &V = cotree.getNode(v);
        if (used[v] == false) {
            used[v] = true;
            for (auto child = V.first_child; child != NetworKit::none;
                    child = cotree.getNode(child).next_sibling) {
                st.push(child);
            }
        } else {
            st.pop();
            if (V.type == NodeType::LEAF) {
                independent_set_size[v] = 1;
            } else if (V.type == NodeType::COMPLEMENT_NODE) {
                for (auto child = V.first_child; child != NetworKit::none;
                        child = cotree.getNode(child).next_sibling) {
                    independent_set_size[v] = std::max(
                        independent_set_size[v], independent_set_size[child]);
                }
            } else {
                for (auto child = V.first_child; child != NetworKit::none;
                        child = cotree.getNode(child).next_sibling) {
                    independent_set_size[v] += independent_set_size[child];
                }
            }
        }
    }
}

void CographIndependentSet::add_to_set() {
    while (!st.empty()) {
        NetworKit::node v = st.top();
        Conode &V = cotree.getNode(v);
        st.pop();
        if (V.type == NodeType::LEAF) {
            independentSet.push_back(v);
        } else if (V.type == NodeType::COMPLEMENT_NODE) {
            NetworKit::node best = NetworKit::none;
            for (auto child = V.first_child; child != NetworKit::none;
                    child = cotree.getNode(child).next_sibling) {
                if (best == NetworKit::none
                        || independent_set_size[child] > independent_set_size[best]) {
                    best = child;
                }
            }
            if (best != NetworKit::none) {
                st.push(best);
            }
        } else {
            for (auto child = V.first_child; child != NetworKit::none;
                    child = cotree.getNode(child).next_sibling) {
                st.push(child);
            }
        }
    }
}

void CographIndependentSet::run() {
    hasRun = true;
    independent_set_size.resize(cotree.upperNodeIdBound(), 0);
    st.push(cotree.getRoot());
    used.resize(cotree.upperNodeIdBound(), false);
    recurse_run();
    st.push(cotree.getRoot());
    used.resize(cotree.upperNodeIdBound(), false);
    add_to_set();
}

NetworKit::count CographIndependentSet::bruteForceIndependetSetSize(NetworKit::Graph &Graph) {
    NetworKit::count n = Graph.numberOfNodes(), ans = 0, flag = 0;
    std::vector<NetworKit::count> st(31), independet_set_nodes;
    st[0] = 1;
    for (int i = 1; i <= 30; i++) {
        st[i] = st[i - 1] * 2;
    }
    for (NetworKit::count mask = 1; mask < st[n]; mask++) {
        independet_set_nodes.clear();
        for (int i = 0; i <= 30; i++) {
            if ((st[i] & mask) > 0) {
                independet_set_nodes.push_back(i);
            }
        }
        flag = 0;
        for (std::size_t i = 0; i < independet_set_nodes.size(); i++) {
            for (std::size_t j = i + 1; j < independet_set_nodes.size(); j++) {
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
