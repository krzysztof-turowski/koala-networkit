/*
 * PerfectGraphRecognition.hpp
 *
 *  Created on: 11.11.2021
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

#include <traversal/BFS.hpp>

namespace Koala {

/**
 * @ingroup recognition
 * The class for recognition of perfect graphs procedure from
 * Chudnovsky, Cornuejols, Liu, Seymour, Vuskovic, "Recognizing Berge graphs".
 *
 */
class PerfectGraphRecognition : public NetworKit::Algorithm {

public:
    /**
     * Given an input graph, set up the perfect graph recognition.
     *
     * @param graph The input graph.
     */
    PerfectGraphRecognition(const NetworKit::Graph &graph);

    /**
     * Execute the perfect graph recognition procedure.
     */
    void run();

    /**
     * Return the result found by the algorithm.
     *
     * @return true if the graph is perfect, false otherwise.
     */
    bool isPerfect() const;

private:
    const std::optional<NetworKit::Graph> graph;
    bool is_perfect;
};

} /* namespace Koala */
