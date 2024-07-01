/*
 * CographRecognition.hpp
 *
 *  Created on: 2024
 *      Author: fixikmila
 */
// Copyright 2024 milana

#pragma once

#include "networkit/base/Algorithm.hpp"
#include "networkit/components/ConnectedComponents.hpp"
#include "networkit/graph/Graph.hpp"
#include "networkit/graph/GraphTools.hpp"

#include <optional>
#include <vector>

#include "recognition/CoTree.hpp"

namespace Koala {

class CographRecognition : public NetworKit::Algorithm {
 public:
    explicit CographRecognition(const NetworKit::Graph &graph);
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
    void check() const;

 protected:
    NetworKit::Graph graph;
};

/**
 * @ingroup recognition
 * The class for recognition of cographs procedure from
 * Bretscher, Corneil, Habib, Paul, "A Simple Linear Time LexBFS Cograph Recognition Algorithm".
 *
 */
class BretscherCorneilHabibPaulCographRecognition : public CographRecognition {
 public:
    using CographRecognition::CographRecognition;
    void run();
    bool isCograph() const;

 private:
    class Info {
     public:
        std::vector<std::vector<std::pair<int, int>>> borders;
        std::vector<NetworKit::node> ans;
    };
    Info info;
    void lex_bfs_minus(bool is_complement, std::vector<NetworKit::node> &a);
    bool neighbourhood_subset_property(bool is_complement, std::vector<NetworKit::node> a,
                                       std::vector<std::vector<std::pair<int, int>>> borders);
    bool is_cograph = false;
};

/**
 * @ingroup recognition
 * The class for recognition of cographs procedure from
 * Corneil, Stewart, Perl, "A linear recognition algorithm for cographs".
 *
 */
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
    void run();
    bool isCograph() const;

 private:
    CorneilStewartPerlCographRecognition::State recognition();
    void unmark();
    void mark(CoNode *x);
    std::pair<CoNode*, CorneilStewartPerlCographRecognition::State>find_lowest() const;
    void insert_x_to_cotree(CoNode *u, CoNode *x);
    State is_cograph = State::UNKNOWN;
    CoTree T;
    int mark_count = 0;
    int mark_and_unmarked_count = 0;
    int mark_ever_count = 0;
    std::queue<CoNode*> marked_with_d_equal_to_md;  // TODO(fixikmila): get rid of this
};

/**
 * @ingroup recognition
 * The class for recognition of cographs procedure from
 * Dahlhaus, "Efficient parallel recognition algorithms of cographs
 * and distance hereditary graphs".
 *
 */
class DahlhausCographRecognition : public CographRecognitionMilana {
 public:
    using CographRecognitionMilana::CographRecognitionMilana;

    void run();

    bool isCograph() const;

 private:
    const int A = 10;
    std::vector<CoNode *> pointer;
    std::vector<CoTree> save;
    bool is_cograph = false;

    CoTree &build_cotree(NetworKit::Graph G, std::vector<int> real_number_of_node);

    void high_low_case(CoTree &T, NetworKit::Graph &G, std::vector<int> &real_number_of_node);
    
    void big_component(CoTree &T, NetworKit::Graph &G, std::vector<int> &vec, std::vector<int> &real_number_of_node);

    inline void add(int vertex_type, CoTree &T, std::vector<int> &vec,
        std::vector<int> &fake_number_of_node, NetworKit::Graph &G,
        std::vector<int> &real_number_of_node);

    bool check_cotree(CoTree T);
};

} /* namespace Koala */
