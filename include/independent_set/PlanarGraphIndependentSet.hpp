/*
 * PlanarGraphIndependentSet.hpp
 *
 *  Created on: 06.04.2023
 *      Author: Mikołaj Twaróg
 *      Ported by: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

#include <independent_set/IndependentSet.hpp>

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
 * The class for the planar graph independent set PTAS using Baker technique.
 *
 */
class BakerPlanarGraphIndependentSet : public PlanarGraphIndependentSet {
 public:
    BakerPlanarGraphIndependentSet(NetworKit::Graph &graph, double epsilon);

    /**
     * Execute the Baker planar graph independent set PTAS.
     */
    void run();

 private:
    double epsilon;
};

/**
 * @ingroup independent_set
 * The class for the planar graph independent set PTAS using Bodlaender technique.
 *
 */
class BodlaenderPlanarGraphIndependentSet : public PlanarGraphIndependentSet {
 public:
    BodlaenderPlanarGraphIndependentSet(NetworKit::Graph &graph, double epsilon);

    /**
     * Execute the Bodlaender planar graph independent set PTAS.
     */
    void run();

 private:
    double epsilon;
};

}  /* namespace Koala */
