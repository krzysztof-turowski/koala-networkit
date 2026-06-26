/*
 * ChordalCliqueCover.cpp
 *
 *  Created on: 2026-06-26
 *      Author: Mateusz Przebieracz
 */

#include <utility>
#include <vector>

#include "clique_cover/ChordalCliqueCover.hpp"

namespace Koala {

ChordalMinCliqueCover::ChordalMinCliqueCover(
    const NetworKit::Graph &graph, const PerfectEliminationOrdering &peo)
    : MinCliqueCover(graph), peo(peo) {
}

void ChordalMinCliqueCover::run() {
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
    clique_cover.clear();

    for (NetworKit::index i = 1; i <= n; i++) {
        NetworKit::node v = peo->alpha_inv[i];

        if (!covered[v]) {
            std::vector<NetworKit::node> C_v;
            C_v.push_back(v);
            covered[v] = true;

            for (auto w : graph->neighborRange(v)) {
                if (peo->alpha[w] > i && !covered[w]) {
                    C_v.push_back(w);
                    covered[w] = true;
                }
            }

            clique_cover.push_back(std::move(C_v));
        }
    }
}

} /* namespace Koala */
