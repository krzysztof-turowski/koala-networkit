/*
 * PerfectGraphVertexColoring.hpp
 *
 *  Created on: 30.03.2023
 *      Author: Adrian Siwiec
 *      Ported by: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

#include <set>

#include "VertexColoring.hpp"

namespace Koala {

/**
 * @ingroup coloring
 * The class for the perfect graph vertex coloring algorithm using SDP.
 *
 */
class PerfectGraphVertexColoring : public VertexColoring {

public:
    using VertexColoring::VertexColoring;

    /**
     * Execute the Groetschel-Lovasz-Schrijver perfect graph vertex coloring algorithm using SDP.
     */
    void run();

    /**
     * Verify the result found by the algorithm.
     */
    void check() const;
private:
    int omega;

    std::vector<int> get_stable_set_intersecting_all_maximum_cliques();

    static int get_theta(const NetworKit::Graph&, const std::vector<int>&);
    static int get_omega(const NetworKit::Graph&);
    static std::vector<int> get_maximum_clique(const NetworKit::Graph&);
    static std::vector<int> get_maximum_stable_set(const NetworKit::Graph&);
    static std::vector<int> get_maximum_weighted_stable_set(
        const NetworKit::Graph&, const std::vector<int>&);
};

} /* namespace Koala */
