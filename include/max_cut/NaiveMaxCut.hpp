/*
 * NaiveMaxCut.hpp
 *
 * Header file for the Max-Cut problem solver using a naive approach.
 * Created on: 13.05.2024
 * Author: Michał Miziołek
 */

#pragma once

#include <vector>
#include <limits>

namespace Koala {

/**
 * @ingroup graph_algorithms
 * The class for solving the Max-Cut problem on a given graph
 * using a naive algorithm.
 */
class NaiveMaxCut {
public:
    /**
     * Constructor that initializes the solver with a graph.
     *
     * @param graphInput The graph on which to solve the Max-Cut problem, represented as an adjacency matrix.
     */
    NaiveMaxCut(const std::vector<std::vector<int>>& graphInput);

    /**
     * Executes the Max-Cut problem solver.
     */
    void solve();

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
    const std::vector<bool>& getBestSet() const;

private:
    std::vector<std::vector<int>> graph;
    int numberOfVertices;
    int maxCutValue;
    std::vector<bool> bestSet;

    // Helper functions
    int calculateCutValue(const std::vector<bool>& set);
};

}  // namespace Koala
