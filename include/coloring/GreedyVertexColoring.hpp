/*
 * GreedyVertexColoring.hpp
 *
 *  Created on: 22.10.2021
 *      Author: Krzysztof Turowski
 */

#ifndef KOALA_COLORING_GREEDY_VERTEX_COLORING_HPP_
#define KOALA_COLORING_GREEDY_VERTEX_COLORING_HPP_

#include <map>

#include <networkit/graph/Graph.hpp>

namespace Koala {

/**
 * @ingroup coloring
 * A set of greedy vertex coloring algorithms.
 *
 */
class GreedyVertexColoring final {

public:
    GreedyVertexColoring() = default;

    /**
     * Given an input graph, find the greedy coloring using Random Sequential algorithm.
     *
     * @param[in]  graph   input graph
     * @param[in]  colors  mapping of colors
     * @param[out]  the number of colors used by the coloring
     */
    int randomSequential(const NetworKit::Graph &graph, std::map<NetworKit::node, int> &colors);

    /**
     * Given an input graph, find the greedy coloring using Largest First algorithm.
     *
     * @param[in]  graph   input graph
     * @param[in]  colors  mapping of colors
     * @param[out]  the number of colors used by the coloring
     */
    int largestFirst(const NetworKit::Graph &graph, std::map<NetworKit::node, int> &colors);

private:
    std::map<NetworKit::node, int>::iterator greedy_color(
        const NetworKit::Graph &graph, NetworKit::node v, std::map<NetworKit::node, int> &colors);
    std::vector<NetworKit::node> largest_first_ordering(
        const NetworKit::Graph &G, std::map<NetworKit::node, int> &colors);
};

} /* namespace Koala */
#endif // KOALA_COLORING_GREEDY_VERTEX_COLORING_HPP_
