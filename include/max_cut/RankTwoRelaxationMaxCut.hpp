/*
 * RankTwoRelaxationMaxCut.hpp
 *
 * Header file for the Max-Cut problem solver using the Rank Two Relaxation method.
 * Created on: 05.04.2024
 * Author: Michał Miziołek
 */

#pragma once

#include <vector>

namespace Koala {

/**
 * @ingroup graph_algorithms
 * The class for solving the Max-Cut problem on a given graph
 * using the Rank Two Relaxation algorithm.
 */
class RankTwoRelaxationMaxCut {
public:
    /**
     * Constructor that initializes the solver with a graph.
     *
     * @param graphInput The graph on which to solve the Max-Cut problem, represented as an adjacency matrix.
     */
    RankTwoRelaxationMaxCut(const std::vector<std::vector<int>>& graphInput);

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
     * 1 represents vertices in one set, while -1 represents vertices in the opposite set.
     *
     * @return Vector of integers indicating the set membership of each vertex.
     */
    const std::vector<int>& getBestSet() const;

private:
    std::vector<std::vector<int>> graph; // Graph represented as an adjacency matrix
    int numberOfVertices;
    double bestCutValue;
    std::vector<int> bestSet;
    std::vector<double> theta;

    // Helper functions
    void distributeThetaEvenly();
    double computeCutValue(const std::vector<int>& x);
    std::vector<int> procedureCut();
    std::vector<double> calculateGradient(const std::vector<double>& theta);
    void gradientDescent(double alpha, int maxIterations);
    void perturbTheta();
};

}  // namespace Koala
