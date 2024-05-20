/*
 * BranchAndBoundMaxCut.hpp
 *
 * Header file for the Max-Cut problem solver using the Branch and Bound method.
 * Created on: 26.03.2024
 * Author: Michał Miziołek
 */

#pragma once

#include <vector>

#include "MaxCut.hpp"

namespace Koala {

/**
 * @ingroup max-cut
 * The class for solving the Max-Cut problem on a given graph
 * using the branch and bound algorithm.
 */
class BranchAndBoundMaxCut final : public MaxCut {
 public:
    using MaxCut::MaxCut;

    /**
     * Executes the Max-Cut problem solver.
     */
    void run();

 private:
    struct Node;

    // Helper functions
    int bound(Node u);
    void branchAndBound();
};

}  // namespace Koala
