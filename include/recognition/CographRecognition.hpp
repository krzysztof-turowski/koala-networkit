/*
 * CographRecognition.hpp
 *
 *  Created on: 29.10.2023
 *      Author: fixikmila
 */

#pragma once

#include <optional>
#include <utility>
#include <vector>

#include <networkit/base/Algorithm.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/graph/AdjListGraph.hpp>
#include <networkit/graph/GraphTools.hpp>

#include "structures/Cotree.hpp"

namespace Koala {

class CographRecognition : public NetworKit::Algorithm {
 public:
    enum class State {
        UNKNOWN,
        COGRAPH,
        NOT_COGRAPH
     };

    /**
     * Given an input graph, set up the cograph recognition.
     *
     * @param graph The input graph.
     */
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
    virtual bool isCograph() const;

    /**
     * Return the graph type found by the algorithm.
     *
     * @return State of the graph.
     */
    virtual State getState() const;

    /**
     * Verify the result found by the algorithm.
     */
    void check() const;

 protected:
    NetworKit::Graph graph;
    State is_cograph;
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

    /**
     * Execute the Bretscher-Corneil-Habib-Paul cograph recognition algorithm.
     */
    void run();
 private:
    class Info {
     public:
        std::vector<std::vector<NetworKit::count>> slice_starts;
        std::vector<NetworKit::node> order;
    };

    Info lex_bfs_minus(bool is_complement, const std::vector<NetworKit::node> &order);
    bool neighbourhood_subset_property(
        bool is_complement, const std::vector<NetworKit::node> &order,
        const std::vector<std::vector<NetworKit::count>> &slice_starts);
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

    /**
     * Execute the Corneil-Stewart-Perl cograph recognition algorithm.
     */
    void run();

 private:
    enum class Marked {
        UNMARKED,
        MARKED,
        MARKED_AND_UNMARKED
    };

    std::vector<NetworKit::node> mark(
        const std::vector<NetworKit::node> &processed_vertices);
    NetworKit::node find_lowest(const std::vector<NetworKit::node> &marked_nodes);
    void attach_to_cotree(NetworKit::node u, NetworKit::node x);

    Cotree T;
    std::vector<Marked> status;
    std::vector<NetworKit::count> md;
    std::vector<NetworKit::node> touched;
};

/**
 * @ingroup recognition
 * The class for recognition of cographs procedure from
 * Dahlhaus, "Efficient parallel recognition algorithms of cographs
 * and distance hereditary graphs".
 *
 */
class DahlhausCographRecognition : public CographRecognition {
 public:
    using CographRecognition::CographRecognition;

    /**
     * Execute the sequential version of Dahlhaus' cograph recognition algorithm.
     */
    void run();

 private:
    const NetworKit::count A = 10;
    Cotree T;
    std::vector<NetworKit::node> covertex;

    NetworKit::node build_cotree(NetworKit::Graph &G);
    void high_low_case(NetworKit::Graph &G);
    void big_component(NetworKit::Graph &G, std::vector<NetworKit::node> &component_nodes);
    inline void attach_to_cotree(
        NodeType node_type, NetworKit::Graph &G, std::vector<NetworKit::node> &subtree_nodes);
    bool check_cotree();
};

class HabibPaulCographRecognition : public CographRecognition {
 public:
    using CographRecognition::CographRecognition;

    void run();

    Cotree cotree;
};

} /* namespace Koala */
