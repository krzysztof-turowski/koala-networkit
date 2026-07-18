/*
 * CliqueCover.hpp
 *
 *  Created on: 2026-06-26
 *      Author: Mateusz Przebieracz
 */

#pragma once

#include <optional>
#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/AdjListGraph.hpp>

namespace Koala {

/**
 * @ingroup clique_cover
 * The base class for the minimum clique cover algorithms.
 */
class MinCliqueCover : public NetworKit::Algorithm {
 public:
    /**
     * Given an input graph, set up the minimum clique cover procedure.
     *
     * @param graph The input graph.
     */
    explicit MinCliqueCover(const NetworKit::Graph &graph);

    /**
     * Return the clique cover found by the algorithm.
     *
     * @return a vector of vectors, where each inner vector represents a clique.
     */
    const std::vector<std::vector<NetworKit::node>>& getCliqueCover() const;

    virtual void run() = 0;

 protected:
    std::optional<NetworKit::Graph> graph;
    std::vector<std::vector<NetworKit::node>> clique_cover;
};

} /* namespace Koala */
