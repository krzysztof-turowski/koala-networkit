/*
 * KargerMinCut.hpp
 *
 * Header file for the Min-Cut problem solver using Karger's algorithm.
 * Created on: 03.05.2024
 * Author: Michał Miziołek
 */

#pragma once

#include <vector>

#include "MinCut.hpp"

namespace Koala {

/**
 * @ingroup min-cut
 * The class for solving the Min-Cut problem on a given graph
 * using Karger's algorithm.
 */
class KargerMinCut final : public MinCut {
 public:
    using MinCut::MinCut;

    /**
     * Executes the Min-Cut problem solver `repeat` times and takes the best solution.
     */
    void run();

    /**
     * Executes the Min-Cut problem solver without repeat.
     */
    void runOnce();

 private:
    int repeat = 10;

    // Helper functions
    int find(std::vector<int>& parent, int i);
    void unionSub(std::vector<int>& parent, std::vector<int>& rank, int x, int y);
};

}  // namespace Koala
