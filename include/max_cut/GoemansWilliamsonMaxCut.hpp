/*
 * GoemansWilliamsonMaxCut.hpp
 *
 * Header file for the Max-Cut problem solver using the Goemans-Williamson algorithm.
 * Created on: 22.04.2024
 * Author: Michał Miziołek
 */

#pragma once

#include <vector>

#include "MaxCut.hpp"

namespace Koala {

/**
 * @ingroup max-cut
 * The class for solving the Max-Cut problem on a given graph
 * using the Goemans-Williamson algorithm.
 */
class GoemansWilliamsonMaxCut final : public MaxCut {
 public:
    using MaxCut::MaxCut;

    /**
     * Executes the Max-Cut problem solver.
     */
    void run();

 private:
    // Helper functions
    std::vector<double> randomUnitVector(int dim);
    void initializeSDP(struct blockmatrix &C, double *&b, struct constraintmatrix *&constraints);
};

}  // namespace Koala
