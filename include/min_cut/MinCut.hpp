/*
 * MinCut.hpp
 *
 * Header file for the MinCut.hpp.
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
 * @ingroup min-cut
 * The base class for the min-cut algorithms.
 *
 */
class MinCut : public NetworKit::Algorithm {
 public:
    /**
     * Given an input graph, set up the min-cut procedure.
     *
     * @param graph The input graph.
     */
    explicit MinCut(NetworKit::Graph &graphInput);

    /**
     * Retrieves the minimum cut value found by the algorithm.
     *
     * @return The minimum cut value.
     */
    int getMinCutValue() const;

    /**
     * Retrieves the best set partition found by the algorithm.
     * True represents vertices in one set, while false represents vertices in the opposite set.
     *
     * @return Vector of boolean indicating the set membership of each vertex.
     */
    const std::vector<bool>& getMinCutSet() const;

 protected:
    std::optional<NetworKit::Graph> graph;
    double minCutValue;
    std::vector<bool> minCutSet;

    // Helper functions
    double calculateCutValue(const std::vector<bool>& set);
};

}  // namespace Koala
