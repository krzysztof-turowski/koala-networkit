/*
 * ChordalClique.cpp
 *
 *  Created on: 2026-06-26
 *      Author: Mateusz Przebieracz
 */

#include "clique/ChordalClique.hpp"

namespace Koala {

ChordalMaxClique::ChordalMaxClique(
    const NetworKit::Graph &graph, const PerfectEliminationOrdering &peo)
    : MaxClique(graph), peo(peo) {
}

void ChordalMaxClique::run() {
    hasRun = true;
    NetworKit::count n = graph->numberOfNodes();

    if (n == 0) return;

    if (!peo.has_value()) {
        MaximumCardinalitySearchChordalGraphRecognition recognizer(*graph);
        recognizer.run();
        peo = recognizer.getPEO();
    }

    NetworKit::count omega_max = 0;
    NetworKit::node v_max = NetworKit::none;

    for (NetworKit::index i = n; i >= 1; i--) {
        NetworKit::node v = peo->alpha_inv[i];
        NetworKit::count n_plus_size = 0;

        for (auto w : graph->neighborRange(v)) {
            if (peo->alpha[w] > i) {
                n_plus_size++;
            }
        }

        if (n_plus_size + 1 > omega_max) {
            omega_max = n_plus_size + 1;
            v_max = v;
        }
    }

    max_clique.clear();
    if (v_max != NetworKit::none) {
        max_clique.push_back(v_max);
        for (auto w : graph->neighborRange(v_max)) {
            if (peo->alpha[w] > peo->alpha[v_max]) {
                max_clique.push_back(w);
            }
        }
    }
}

} /* namespace Koala */
