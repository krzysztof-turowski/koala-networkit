/*
 * BranchAndBoundMaxCut.cpp
 *
 * Solution for the Max-Cut problem using Branch and Bound method.
 * Created on: 26.03.2024
 * Author: Michał Miziołek
 */

#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>

#include <max_cut/BranchAndBoundMaxCut.hpp>

namespace Koala {

struct BranchAndBoundMaxCut::Node {
    std::vector<bool> set;      // Indicates which set a vertex belongs to
    int level;                  // Current vertex to consider
    int bound;                  // Upper bound of the maximum cut value from this node
};

int BranchAndBoundMaxCut::bound(Node u) {
    int result = calculateCutValue(u.set);
    graph->forEdges([&](NetworKit::node j, NetworKit::node k, NetworKit::edgeweight w) {
        if (j >= u.level) {
            result += w;
        }
    });
    return result;
}

void BranchAndBoundMaxCut::branchAndBound() {
    maxCutValue = 0;
    std::vector<Node> stack;
    Node root;
    root.level = 0;
    root.set.resize(graph->numberOfNodes(), false);
    root.bound = bound(root);
    stack.push_back(root);

    while (!stack.empty()) {
        Node u = stack.back();
        stack.pop_back();

        if (u.level == graph->numberOfNodes()) {
            double currentCutValue = calculateCutValue(u.set);
            if (currentCutValue > maxCutValue) {
                maxCutValue = currentCutValue;
                maxCutSet = u.set;
            }
        } else {
            for (int i = 0; i < 2; ++i) {
                Node v = u;
                v.level = u.level + 1;
                v.set[u.level] = i;
                v.bound = bound(v);
                if (v.bound > maxCutValue) {
                    stack.push_back(v);
                }
            }
        }
    }
}

void BranchAndBoundMaxCut::run() {
    branchAndBound();
}

}  // namespace Koala
