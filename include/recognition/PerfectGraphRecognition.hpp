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
    enum class State {
        UNKNOWN,
        PERFECT,
        HAS_JEWEL,
        HAS_PYRAMID,
        HAS_T1,
        HAS_T2,
        HAS_T3,
        HAS_NEAR_CLEANER_ODD_HOLE
     };

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

    /**
     * Return the graph type found by the algorithm.
     *
     * @return State of the graph.
     */
    State getState() const;

private:
    const std::optional<NetworKit::Graph> graph;
    State is_perfect;

    static State containsSimpleProhibited(const NetworKit::Graph &graph);
    static bool containsJewel(const NetworKit::Graph &graph);
};

bool PerfectGraphRecognitionNaive(const NetworKit::Graph &graph);

} /* namespace Koala */
