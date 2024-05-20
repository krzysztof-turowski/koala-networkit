/*
 * HaoOrlinMinCut.hpp
 *
 * Header file for the Min-Cut problem solver using the Hao-Orlin algorithm.
 * Created on: 30.04.2024
 * Author: Michał Miziołek
 */

#pragma once

#include <vector>

#include "MinCut.hpp"

namespace Koala {

/**
 * @ingroup min-cut
 * The class for solving the Min-Cut problem on a given graph
 * using the Hao-Orlin algorithm.
 */
template <typename MaxFlowT>
class HaoOrlinMinCut final : public MinCut {
 public:
    using MinCut::MinCut;

    /**
     * Executes the Min-Cut problem solver `repeat` times and takes the best solution.
     */
    void run() {
        minCutValue = INT_MAX;
        minCutSet.assign(graph->numberOfNodes(), false);
        std::vector<bool> visited(graph->numberOfNodes(), false);
        visited[0] = true;
        std::vector<int> S = {0};

        while (S.size() < graph->numberOfNodes()) {
            int t_prime = -1;
            for (int i = 0; i < graph->numberOfNodes(); ++i) {
                if (!visited[i]) {
                    t_prime = i;
                    break;
                }
            }

            MaxFlowT minCutSolver(*graph, 0, t_prime);
            minCutSolver.run();
            double currentMinCutValue = minCutSolver.getFlowSize();

            minCutValue = std::min(minCutValue, currentMinCutValue);

            visited[t_prime] = true;
            S.push_back(t_prime);
        }
    }

 private:
    // Helper functions
};

}  // namespace Koala
