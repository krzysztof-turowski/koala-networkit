/*
 * PerfectGraphRecognition.hpp
 *
 *  Created on: 11.11.2021
 *      Author: Adrian Siwiec
 *      Ported by: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

#include <optional>
#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/GraphTools.hpp>

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
    explicit PerfectGraphRecognition(NetworKit::Graph &graph);

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

    /**
     * Verify the result found by the algorithm.
     */
    void check() const;

    static bool isComplete(
        const NetworKit::Graph &graph, const std::vector<NetworKit::node> &X, NetworKit::node v);
    static std::vector<NetworKit::node> getAllCompleteVertices(
        const NetworKit::Graph &graph, const std::vector<NetworKit::node> &X);
    static std::vector<std::vector<NetworKit::node>> getAuxiliaryComponents(
        const NetworKit::Graph &graph, const std::vector<NetworKit::node> &V);
 private:
    std::optional<NetworKit::Graph> graph;
    State is_perfect;

    static State contains_simple_prohibited(const NetworKit::Graph &graph);
    static bool contains_jewel(const NetworKit::Graph &graph);
    static bool contains_pyramid(const NetworKit::Graph &graph);
    static bool contains_t1(const NetworKit::Graph &graph);
    static bool contains_t2(const NetworKit::Graph &graph);
    static bool contains_t3(const NetworKit::Graph &graph);
    static bool contains_near_cleaner_odd_hole(const NetworKit::Graph &graph);

    static bool contains_odd_hole(const NetworKit::Graph &graph);
    static bool contains_hole(const NetworKit::Graph &graph, NetworKit::count length);
};

} /* namespace Koala */
