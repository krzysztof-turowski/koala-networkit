/*
 * PlanarGraphRecognition.hpp
 *
 *  Created on: 06.04.2023
 */

#pragma once

#include <optional>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace Koala {

/**
 * @ingroup recognition
 *
 */
class PlanarGraphRecognition : public NetworKit::Algorithm {
 public:
    /**
     * Given an input graph, set up the planar graph recognition.
     *
     * @param graph The input graph.
     */
    PlanarGraphRecognition(NetworKit::Graph &graph);

    /**
     * Execute the planar graph recognition procedure.
     */
    void run();

    /**
     * Return the result found by the algorithm.
     *
     * @return true if the graph is planar, false otherwise.
     */
    bool isPlanar() const;

 private:
    std::optional<NetworKit::Graph> graph;
    bool is_planar;
};

}  /* namespace Koala */
