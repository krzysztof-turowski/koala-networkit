/*
 * StoerWagnerMinCut.hpp
 *
 * Header file for the Min-Cut problem solver using the Stoer-Wagner algorithm.
 * Created on: 07.05.2024
 * Author: Michał Miziołek
 */

#pragma once

#include <vector>

#include "MinCut.hpp"

namespace Koala {

/**
 * @ingroup min-cut
 * The class for solving the Min-Cut problem on a given graph
 * using the Stoer-Wagner algorithm.
 */
class StoerWagnerMinCut final : public MinCut {
 public:
    using MinCut::MinCut;

    /**
     * Executes the Min-Cut problem solver.
     */
    void run();

 private:
    // Helper functions
    double minCutPhase(std::vector<NetworKit::node>& vertices,
                        std::vector<std::vector<double>>& edges);
};

}  // namespace Koala
