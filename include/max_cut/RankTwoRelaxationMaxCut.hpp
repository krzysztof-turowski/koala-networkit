/*
 * RankTwoRelaxationMaxCut.hpp
 *
 * Header file for the Max-Cut problem solver using the Rank Two Relaxation method.
 * Created on: 05.04.2024
 * Author: Michał Miziołek
 */

#pragma once

#include <vector>

#include "MaxCut.hpp"

namespace Koala {

/**
 * @ingroup max-cut
 * The class for solving the Max-Cut problem on a given graph
 * using the Rank Two Relaxation algorithm.
 */
class RankTwoRelaxationMaxCut final : public MaxCut {
 public:
    using MaxCut::MaxCut;

    /**
     * Executes the Max-Cut problem solver.
     */
    void run();

 private:
    static const double alpha = 0.001;
    static const int maxIterations = 100000;
    std::vector<double> theta;

    // Helper functions
    void distributeThetaEvenly();
    std::vector<bool> procedureCut();
    std::vector<double> calculateGradient(const std::vector<double>& theta);
    void gradientDescent(double alpha, int maxIterations);
    void perturbTheta();
};

}  // namespace Koala
