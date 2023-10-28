

#pragma once

#include <optional>
#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace Koala {


class CoraphRecognition : public NetworKit::Algorithm {
 public:
    
    enum class State {
            UNKNOWN,
            COMPLEMENT_REDUCIBLE,
            COND_1,
            COND_2,
            COND_3,
            COND_4,
            COND_5,
            COND_6
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
     * @return true if the graph is complement reducible, false otherwise.
     */
    bool isComplementReducible() const;
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

    static bool isComplementReducible(
        const NetworKit::Graph &graph, const std::vector<NetworKit::node> &X, NetworKit::node v);
 private:
    std::optional<NetworKit::Graph> graph;
    State is_complement_reducible;


};

} /* namespace Koala */
