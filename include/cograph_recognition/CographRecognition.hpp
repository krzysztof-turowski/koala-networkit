#pragma once

#include <optional>
#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace Koala {


class CographRecognition : public NetworKit::Algorithm {
 public:
    
    enum class State {
            UNKNOWN,
            COGRAPH,
            CONTAINS_0_NODE,
            EXISTS_1_NODE_NOT_PROPERLY_MARKED,
            GRANDPARENT_IS_NOT_IN_SET,
            NO_ONE_PATH,
            WRONG_PARENT,
            WRONG_GRANDPARENT
         };

    /**
     * Given an input graph, set up the cograph recognition.
     *
     * @param graph The input graph.
     */
    explicit CographRecognition(NetworKit::Graph &graph);

    /**
     * Execute the cograph recognition procedure.
     */
    void run();

    /**
     * Return the result found by the algorithm.
     *
     * @return true if the graph is a cograph, false otherwise.
     */
    bool isCograph() const;
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
    static State Cograph_Recognition(NetworKit::Graph &graph);
 private:
    NetworKit::Graph graph;
    State is_cograph;


};

} /* namespace Koala */
