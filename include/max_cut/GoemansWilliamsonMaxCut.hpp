/*
 * GoemansWilliamsonMaxCut.hpp
 *
 * Header file for the Max-Cut problem solver using the Goemans-Williamson algorithm.
 * Created on: 22.04.2024
 * Author: Michał Miziołek
 */

#pragma once

#include <vector>
#include <iostream>

extern "C" {
#include <declarations.h>
}

namespace Koala {

/**
 * @ingroup graph_algorithms
 * The class for solving the Max-Cut problem on a given graph
 * using the Goemans-Williamson algorithm.
 */
class GoemansWilliamsonMaxCut {
public:
    /**
     * Constructor that initializes the solver with a graph.
     *
     * @param graphInput The graph on which to solve the Max-Cut problem, represented as an adjacency matrix.
     */
    GoemansWilliamsonMaxCut(const std::vector<std::vector<int>>& graphInput);

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
    int n;
    std::vector<bool> bestSet;
    int maxCutValue;

    std::vector<double> randomUnitVector(int dim);
    void initializeSDP(struct blockmatrix &C, double *&b, struct constraintmatrix *&constraints);
    void freeSDP(struct blockmatrix &C, struct blockmatrix &X, struct blockmatrix &Z, double *b, struct constraintmatrix *constraints, double *y);
};

} // namespace Koala
