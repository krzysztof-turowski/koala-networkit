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

void NaiveMaxCut::run() {
    maxCutValue = 0;
    std::vector<bool> set(graph->numberOfNodes(), false);
    bool improved = true;

    while (improved) {
        improved = false;
        for (int i = 0; i < graph->numberOfNodes(); i++) {
            set[i] = !set[i];
            double newCut = calculateCutValue(set);
            if (newCut > maxCutValue) {
                maxCutValue = newCut;
                maxCutSet = set;
                improved = true;
            } else {
                set[i] = !set[i];
            }
        }
    }
}

}  // namespace Koala
