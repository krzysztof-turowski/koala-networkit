/*
 * ChordalCliqueCover.hpp
 *
 *  Created on: 2026-06-26
 *      Author: Mateusz Przebieracz
 */

#pragma once

#include <optional>

#include "clique_cover/CliqueCover.hpp"
#include "recognition/ChordalGraphRecognition.hpp"

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
