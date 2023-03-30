/*
 * GreedyVertexColoring.hpp
 *
 *  Created on: 22.10.2021
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

#include "VertexColoring.hpp"

namespace Koala {

/**
 * @ingroup coloring
 * The base class for the greedy vertex coloring heuristics.
 *
 */
class GreedyVertexColoring : public VertexColoring {

public:
    using VertexColoring::VertexColoring;

protected:
    std::map<NetworKit::node, int>::iterator greedy_color(NetworKit::node v);
};

/**
 * @ingroup coloring
 * The class for the random sequential greedy vertex coloring heuristic.
 */
class RandomSequentialVertexColoring final : public GreedyVertexColoring {

public:
    using GreedyVertexColoring::GreedyVertexColoring;

    /**
     * Execute the random sequential greedy vertex coloring heuristic.
     */
    void run();
};

/**
 * @ingroup coloring
 * The class for the largest first greedy vertex coloring heuristic.
 */
class LargestFirstVertexColoring final : public GreedyVertexColoring {

public:
    using GreedyVertexColoring::GreedyVertexColoring;

    /**
     * Execute the largest first greedy vertex coloring heuristic.
     */
    void run();

private:
    std::vector<NetworKit::node> largest_first_ordering();
};

/**
 * @ingroup coloring
 * The class for the smallest last greedy vertex coloring heuristic.
 */
class SmallestLastVertexColoring final : public GreedyVertexColoring {

public:
    using GreedyVertexColoring::GreedyVertexColoring;

    /**
     * Execute the smallest last greedy vertex coloring heuristic.
     */
    void run();

private:
    std::vector<NetworKit::node> smallest_last_ordering();
};

/**
 * @ingroup coloring
 * The class for the saturated largest first greedy vertex coloring heuristic.
 */
class SaturatedLargestFirstVertexColoring final : public GreedyVertexColoring {

public:
    using GreedyVertexColoring::GreedyVertexColoring;

    /**
     * Execute the saturated largest first greedy vertex coloring heuristic.
     */
    void run();
};

/**
 * @ingroup coloring
 * The class for the greedy independent set vertex coloring heuristic.
 */
class GreedyIndependentSetVertexColoring final : public GreedyVertexColoring {

public:
    using GreedyVertexColoring::GreedyVertexColoring;

    /**
     * Execute the greedy independent set vertex coloring heuristic.
     */
    void run();
};

} /* namespace Koala */
