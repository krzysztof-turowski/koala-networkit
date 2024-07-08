/*
 * MaxCut.hpp
 *
 * Header file for the Max-Cut.hpp.
 * Created on: 17.05.2024
 * Author: Michał Miziołek
 */

#pragma once

#include <map>
#include <optional>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace Koala {

/**
 * @ingroup max-cut
 * The base class for the max-cut algorithms.
 *
 */
class MaxCut : public NetworKit::Algorithm {
 public:
    /**
     * Given an input graph, set up the max-cut procedure.
     *
     * @param graph The input graph.
     */
    explicit MaxCut(NetworKit::Graph &graphInput);

    /**
     * Retrieves the maximum cut value found by the algorithm.
     *
     * @return The maximum cut value.
     */
    int getMaxCutValue() const;

    /**
     * Retrieves the best set partition found by the algorithm.
     * True represents vertices in one set, while false represents vertices in the opposite set.
     *
     * @return Vector of boolean indicating the set membership of each vertex.
     */
    const std::vector<bool>& getMaxCutSet() const;

 protected:
    std::optional<NetworKit::Graph> graph;
    double maxCutValue;
    std::vector<bool> maxCutSet;

    // Helper functions
    double calculateCutValue(const std::vector<bool>& set);
};

}  // namespace Koala
