/*
 * PerfectGraphColoring.hpp
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
class PerfectGraphColoring : public VertexColoring {

public:
    using VertexColoring::VertexColoring;

    /**
     * Execute the perfect graph vertex coloring algorithm using SDP.
     */
    void run();

    int get_omega();
    int get_chi();
private:
    std::vector<NetworKit::node> get_stable_set_intersecting_all_maximum_cliques();
    std::vector<int> get_stable_set_intersecting_maximum_cliques(const std::vector<int>&);
    std::vector<int> get_stable_set_intersecting_maximum_cliques_2(const std::vector<int>&);
};

} /* namespace Koala */
