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
    


    /**
     * Given an input graph, set up the cograph recognition.
     *
     * @param graph The input graph.
     */
    explicit CographRecognition(NetworKit::Graph &graph);

    /**
     * Execute the cograph recognition procedure.
     */
    //virtual void run();

    /**
     * Return the result found by the algorithm.
     *
     * @return true if the graph is a cograph, false otherwise.
     */
    //virtual bool isCograph() const;

    /**
     * Verify the result found by the algorithm.
     */
    //virtual void check() const;

protected:
    NetworKit::Graph graph;


};

} /* namespace Koala */
