/*
 * PerfectGraphColoring.hpp
 *
 *  Created on: 30.03.2023
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
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

private:
    std::set<NetworKit::node> get_stable_set_intersecting_all_maximum_cliques();
};

} /* namespace Koala */