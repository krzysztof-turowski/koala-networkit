/*
 * BranchAndBoundMaxCut.hpp
 *
 * Header file for the Max-Cut problem solver using the Branch and Bound method.
 * Created on: 26.03.2024
 * Author: Michał Miziołek
 */

#pragma once

#include <vector>

namespace Koala {

/**
 * @ingroup graph_algorithms
 * The class for solving the Max-Cut problem on a given graph
 * using the branch and bound algorithm.
 */
class BranchAndBoundMaxCut {
public:
    /**
     * Constructor that initializes the solver with a graph.
     *
     * @param graphInput The graph on which to solve the Max-Cut problem, represented as an adjacency matrix.
     */
    BranchAndBoundMaxCut(const std::vector<std::vector<int>>& graphInput);

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
    struct Node;
    std::vector<std::vector<int>> graph; // Graph represented as an adjacency matrix
    int numberOfVertices;
    int maxCutValue;
    std::vector<bool> bestSet;

    // Helper functions
    int calculateCutValue(const std::vector<bool>& set);
    int bound(Node u);
    void branchAndBound();
};

}  // namespace Koala
