/*
 * PlanarGraphIndependentSet.hpp
 *
 *  Created on: 06.04.2023
 *      Author: Mikołaj Twaróg
 *      Ported by: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

#include "IndependentSet.hpp"

namespace Koala {

/**
 * @ingroup independent_set
 * The class for the planar graph independent set algorithm.
 *
 */
class PlanarGraphIndependentSet : public IndependentSet {
    using IndependentSet::IndependentSet;
};

/**
 * @ingroup independent_set
 * The class for the perfect graph independent set algorithm using Baker technique.
 *
 */
class BakerPlanarGraphIndependentSet : public PlanarGraphIndependentSet {
 public:
    BakerPlanarGraphIndependentSet(NetworKit::Graph &graph, double epsilon);

    /**
     * Execute the Baker planar graph approximation independent set PTAS.
     */
    void run();

 private:
    double epsilon;
};

}  /* namespace Koala */
