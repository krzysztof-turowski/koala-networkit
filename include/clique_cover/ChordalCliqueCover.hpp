/*
 * ChordalMinCliqueCover.hpp
 *
 *  Created on: 2026-05-25
 *      Author: Mateusz Przebieracz
 */

#pragma once

#include "clique_cover/CliqueCover.hpp"
#include "recognition/ChordalGraphRecognition.hpp"

#include <optional>

namespace Koala {

/**
 * @ingroup clique_cover
 * Minimum Clique Cover algorithm for Chordal Graphs.
 */
class ChordalMinCliqueCover final : public MinCliqueCover {
 public:
    using MinCliqueCover::MinCliqueCover;

    /**
     * Constructor injecting a pre-calculated PEO.
     * @param graph The input chordal graph.
     * @param peo A valid perfect elimination ordering.
     */
    ChordalMinCliqueCover(const NetworKit::Graph &graph, const PerfectEliminationOrdering &peo);

    void run() override;

 private:
    std::optional<PerfectEliminationOrdering> peo;
};

} /* namespace Koala */