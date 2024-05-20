/*
 * NaiveMaxCut.hpp
 *
 * Header file for the Max-Cut problem solver using a naive approach.
 * Created on: 13.05.2024
 * Author: Michał Miziołek
 */

#pragma once

#include <vector>

#include "MaxCut.hpp"

namespace Koala {

/**
 * @ingroup max-cut
 * The class for solving the Max-Cut problem on a given graph
 * using a naive algorithm.
 */
class NaiveMaxCut final : public MaxCut {
 public:
    using MaxCut::MaxCut;

    /**
     * Executes the Max-Cut problem solver.
     */
    void run();

 private:
    // Helper functions
};

}  // namespace Koala
