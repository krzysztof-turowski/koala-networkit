/*
 * ChordalVertexColoring.hpp
 *
 *  Created on: 2026-06-26
 *      Author: Mateusz Przebieracz
 */

#pragma once

#include <optional>

#include "coloring/VertexColoring.hpp"
#include "recognition/ChordalGraphRecognition.hpp"

namespace Koala {

/**
 * @ingroup coloring
 * Linear-time Minimum Vertex Coloring algorithm for Chordal Graphs.
 */
class ChordalVertexColoring final : public VertexColoring {
 public:
    using VertexColoring::VertexColoring;

    /**
     * Constructor injecting a pre-calculated PEO.
     * @param graph The input chordal graph.
     * @param peo A valid perfect elimination ordering.
     */
    ChordalVertexColoring(const NetworKit::Graph &graph, const PerfectEliminationOrdering &peo);

    void run() override;

 private:
    std::optional<PerfectEliminationOrdering> peo;
};

} /* namespace Koala */
