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
    std::vector<bool> set; // Indicates which set a vertex belongs to
    int level; // Current vertex to consider
    int bound; // Upper bound of the maximum cut value from this node
};

int BranchAndBoundMaxCut::calculateCutValue(const std::vector<bool>& set) {
    int cutValue = 0;
    for (int i = 0; i < numberOfVertices; ++i) {
        for (int j = i + 1; j < numberOfVertices; ++j) {
            if (set[i] != set[j] && graph[i][j] > 0) {
                cutValue += graph[i][j];
            }
        }
    }
    return cutValue;
}

int BranchAndBoundMaxCut::bound(Node u) {
    int result = calculateCutValue(u.set);
    for (int j = u.level; j < numberOfVertices; ++j) {
        for (int k = 0; k < numberOfVertices; ++k) {
            if (graph[j][k] > 0) {
                result += graph[j][k];
            }
        }
    }
    return result;
}

void BranchAndBoundMaxCut::branchAndBound() {
    std::vector<Node> stack;
    Node root;
    root.level = 0;
    root.set.resize(numberOfVertices, false);
    root.bound = bound(root);
    stack.push_back(root);

    while (!stack.empty()) {
        Node u = stack.back();
        stack.pop_back();

        if (u.level == numberOfVertices) {
            int currentCutValue = calculateCutValue(u.set);
            if (currentCutValue > maxCutValue) {
                maxCutValue = currentCutValue;
                bestSet = u.set;
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

BranchAndBoundMaxCut::BranchAndBoundMaxCut(const std::vector<std::vector<int>>& graphInput)
    : graph(graphInput), numberOfVertices(graphInput.size()), maxCutValue(std::numeric_limits<int>::min()) {}

void BranchAndBoundMaxCut::solve() {
    branchAndBound();
}

int BranchAndBoundMaxCut::getMaxCutValue() const {
    return maxCutValue;
}

const std::vector<bool>& BranchAndBoundMaxCut::getBestSet() const {
    return bestSet;
}

}  // namespace Koala
