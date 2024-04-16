/*
 * IPlanarGraphRecognition.hpp
 *
 *  Created on: 24.03.2024
 *      Author: Dzianis Lahunou
 *      Ported by: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

#include <optional>
#include <vector>
#include <limits.h>
#include <stack>
#include <list>
#include <networkit/base/Algorithm.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace Koala {

/**
 * @ingroup recognition
 * The Interface for class for recognition of planar graphs p.
 *
 */
class IPlanarGraphRecognition : public NetworKit::Algorithm {
 public:
    enum class State {
        NOT_PLANAR,
        PLANAR
     };    
     
      /**
         * Given an input graph, set up the planar graph recognition.
         *
         * @param graph The input graph.
        */
     IPlanarGraphRecognition::IPlanarGraphRecognition(NetworKit::Graph &graph, bool embedding);

    /**
     * Execute the planar graph recognition procedure.
     */
    virtual void run() = 0;

    /**
     * Return the result found by the algorithm.
     *
     * @return true if the graph is planar, false otherwise.
     */
   bool isPlanar() const;


    /**
     * Return the embbending found by the algorithm.
     *
     * @return pointer to embedding if the graph is planar, nullptr otherwise.
     */
   void getEmbedding() const;


 protected:
    NetworKit::Graph graph;
    State is_planar;
    bool embedding;

};

} /* namespace Koala */
