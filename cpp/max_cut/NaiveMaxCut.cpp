/*
 * NaiveMaxCut.cpp
 *
 * Solution for the Max-Cut problem using a naive method.
 * Created on: 13.05.2024
 * Author: Michał Miziołek
 */

#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>

#include <max_cut/NaiveMaxCut.hpp>

namespace Koala {

NaiveMaxCut::NaiveMaxCut(const std::vector<std::vector<int>>& graphInput)
    : graph(graphInput), numberOfVertices(graphInput.size()), maxCutValue(0) {
    bestSet.resize(numberOfVertices, false);
}

int NaiveMaxCut::calculateCutValue(const std::vector<bool>& set) {
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

void NaiveMaxCut::solve() {
    std::vector<bool> set(numberOfVertices, false);
    bool improved = true;
    while (improved) {
        improved = false;
        for (int i = 0; i < numberOfVertices; i++) {
            set[i] = !set[i];
            int newCut = calculateCutValue(set);
            if (newCut > maxCutValue) {
                maxCutValue = newCut;
                bestSet = set;
                improved = true;
            } else {
                set[i] = !set[i];
            }
        }
    }
}

int NaiveMaxCut::getMaxCutValue() const {
    return maxCutValue;
}

const std::vector<bool>& NaiveMaxCut::getBestSet() const {
    return bestSet;
}

}  // namespace Koala
