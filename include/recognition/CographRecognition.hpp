#pragma once

#include <optional>
#include <vector>

#include "networkit/base/Algorithm.hpp"
#include "networkit/components/ConnectedComponents.hpp"
#include "networkit/graph/Graph.hpp"
#include "networkit/graph/GraphTools.hpp"

namespace Koala {


class CographRecognition : public NetworKit::Algorithm {
 public:

    /**
     * Execute the cograph recognition procedure.
     */
    virtual void run() = 0;

    /**
     * Return the result found by the algorithm.
     *
     * @return true if the graph is a cograph, false otherwise.
     */
    virtual bool isCograph() const = 0;

    /**
     * Verify the result found by the algorithm.
     */
    virtual void check() const = 0;

 protected:
    NetworKit::Graph graph;
};

} /* namespace Koala */