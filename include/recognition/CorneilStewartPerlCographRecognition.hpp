/*
 * CorneilStewartPerlCographRecognition.hpp
 *
 *  Created on: 2024
 *      Author: fixikmila
 */
// Copyright 2024 milana

#pragma once

#include <optional>
#include <vector>

#include "networkit/base/Algorithm.hpp"
#include "networkit/components/ConnectedComponents.hpp"
#include "networkit/graph/Graph.hpp"
#include "networkit/graph/GraphTools.hpp"
#include "recognition/CographRecognition.hpp"
#include "recognition/CoTree.hpp"

namespace Koala {
class CorneilStewartPerlCographRecognition : public CographRecognition {
 public:
        using CographRecognition::CographRecognition;
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
        CorneilStewartPerlCographRecognition::State Cograph_Recognition();
        void run() override;
        bool isCograph() const override;
        void Unmark();
        void Mark(CoNode *x);
        std::pair<CoNode*, CorneilStewartPerlCographRecognition::State>FindLowest();
        void InsertXToCoTree(CoNode *u, CoNode *x);

 private:
        State is_cograph = State::UNKNOWN;
        CoTree T;
        NetworKit::Graph G;
        int mark_count = 0;
        int mark_and_unmarked_count = 0;
        int mark_ever_count = 0;
        std::queue<CoNode*> marked_with_d_equal_to_md;  // TODO: get rid of this
};
} /* namespace Koala */
