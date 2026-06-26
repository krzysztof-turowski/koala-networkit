/*
 * ChordalIndependentSet.cpp
 *
 *  Created on: 2026-06-26
 *      Author: Mateusz Przebieracz
 */

#include "independent_set/ChordalIndependentSet.hpp"

#include <vector>

namespace Koala {

ChordalIndependentSet::ChordalIndependentSet(
    const NetworKit::Graph &graph, const PerfectEliminationOrdering &peo)
    : IndependentSet(graph), peo(peo) {
}

void ChordalIndependentSet::run() {
    hasRun = true;
    NetworKit::count n = graph->numberOfNodes();
    NetworKit::count bound = graph->upperNodeIdBound();

    if (n == 0) return;

    if (!peo.has_value()) {
        MaximumCardinalitySearchChordalGraphRecognition recognizer(*graph);
        recognizer.run();
        peo = recognizer.getPEO();
    }

    std::vector<uint8_t> covered(bound, false);
    independentSet.clear();

    for (NetworKit::index i = 1; i <= n; i++) {
        NetworKit::node v = peo->alpha_inv[i];

        if (!covered[v]) {
            independentSet.push_back(v);
            covered[v] = true;

            for (auto w : graph->neighborRange(v)) {
                if (peo->alpha[w] > i) {
                    covered[w] = true;
                }
            }
        }
    }
}

} /* namespace Koala */
