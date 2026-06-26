/*
 * ChordalIndependentSet.hpp
 *
 *  Created on: 2026-06-26
 *      Author: Mateusz Przebieracz
 */

#pragma once

#include <optional>

#include "independent_set/IndependentSet.hpp"
#include "recognition/ChordalGraphRecognition.hpp"

namespace Koala {

/**
 * @ingroup independent_set
 * Linear-time Maximum Independent Set algorithm for Chordal Graphs.
 */
class ChordalIndependentSet final : public IndependentSet {
 public:
    using IndependentSet::IndependentSet;

    /**
     * Constructor injecting a pre-calculated PEO.
     * @param graph The input chordal graph.
     * @param peo A valid perfect elimination ordering.
     */
    ChordalIndependentSet(const NetworKit::Graph &graph, const PerfectEliminationOrdering &peo);

    void run() override;

 private:
    std::optional<PerfectEliminationOrdering> peo;
};

} /* namespace Koala */
