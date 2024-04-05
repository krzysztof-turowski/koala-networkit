#pragma once

#include <optional>
#include <vector>

#include "networkit/base/Algorithm.hpp"
#include "networkit/components/ConnectedComponents.hpp"
#include "networkit/graph/Graph.hpp"
#include "networkit/graph/GraphTools.hpp"
#include "recognition/CographRecognition.hpp"

namespace Koala {
    class CorneilStewartPerlCographRecognition : public CographRecognition {
    public:
        explicit CorneilStewartPerlCographRecognition(NetworKit::Graph &graph);
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
         * Return the graph type found by the algorithm.
         *
         * @return State of the graph.
         */
        State getState() const;

        static State Cograph_Recognition(NetworKit::Graph &graph);
        void run() override;
        bool isCograph() const override;
        void check() const override;
    private:
        NetworKit::Graph graph;
        State is_cograph;
    };
} /* namespace Koala */
