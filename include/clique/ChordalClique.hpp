/*
 * ChordalMaxClique.hpp
 *
 *  Created on: 2026-05-25
 *      Author: Mateusz Przebieracz
 */

#pragma once

#include "clique/Clique.hpp"
#include "recognition/ChordalGraphRecognition.hpp"

#include <optional>

namespace Koala {

/**
 * @ingroup clique
 * Maximum Clique algorithm for Chordal Graphs.
 */
class ChordalMaxClique final : public MaxClique {
 public:
    using MaxClique::MaxClique;

    /**
     * Constructor injecting a pre-calculated PEO.
     * @param graph The input chordal graph.
     * @param peo A valid perfect elimination ordering.
     */
    ChordalMaxClique(const NetworKit::Graph &graph, const PerfectEliminationOrdering &peo);

    void run() override;

 private:
    std::optional<PerfectEliminationOrdering> peo;
};

} /* namespace Koala */