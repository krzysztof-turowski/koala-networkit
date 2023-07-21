/*
 * ExactDominatingSet.hpp
 *
 *  Created on: 01.07.2023
 *      Author: Piotr Kubaty
 */

#pragma once

#include <set>

#include <dominating_set/DominatingSet.hpp>

namespace Koala {

/**
 * @ingroup dominating_set
 * The base class for the dominating set algorithms via set covering.
 *
 */
template<typename SetCoverAlgorithm>
class BranchAndReduceDominatingSet : public DominatingSet {
 public:
    using DominatingSet::DominatingSet;

    /**
     * Execute the exact dominating set procedure via set covering.
     */
    void run() {
        hasRun = true;
        std::vector<std::set<NetworKit::node>> family;
        graph->forNodes([&family, this](NetworKit::node u) {
            std::set<NetworKit::node> neighborhood;
            neighborhood.insert(u);
            graph->forNeighborsOf(u, [&neighborhood](NetworKit::node v) {
                neighborhood.insert(v);
            });
            family.emplace_back(neighborhood);
        });
        std::vector<std::set<NetworKit::index>> occurences(family);
        auto set_cover_algorithm = SetCoverAlgorithm(family, occurences);
        set_cover_algorithm.run();
        dominating_set = set_cover_algorithm.getSetCover();
    }
};

/**
 * @ingroup dominating_set
 * The base class for the exact dominating set algorithms.
 *
 */
class ExactDominatingSet : public DominatingSet {
 public:
    using DominatingSet::DominatingSet;

 protected:
    std::set<NetworKit::node> free, bounded, required;
};

/**
 * @ingroup dominating_set
 * The class for the Fomin-Kratsch-Woeginger exact dominating set algorithm.
 *
 */
class FominKratschWoegingerDominatingSet : public ExactDominatingSet {
 public:
    using ExactDominatingSet::ExactDominatingSet;

    /**
     * Execute the Fomin-Kratsch-Woeginger exact dominating set algorithm.
     */
    void run();

 protected:
    std::set<NetworKit::node> degreeOne, degreeTwo;
    std::vector<std::set<NetworKit::node>> neighborhood;

    std::vector<bool> recurse();
    ////////////
    void onDegreeDecrement(NetworKit::node updated);
    void onDegreeIncrement(NetworKit::node updated);
    void forgetVertexInNeighbors(NetworKit::node vertex);
    void retrieveVertexInNeighbors(NetworKit::node vertex);
    bool forgetVertex(NetworKit::node forced);
    void retrieveVertex(NetworKit::node vertex, bool isFree);
    std::tuple<std::vector<NetworKit::node>, bool> addToTheSolution(NetworKit::node forced);
    void removeFromTheSolution(
        NetworKit::node vertex,
        std::vector<NetworKit::node> &movedVertices,
        bool isVertexFree);
};

/**
 * @ingroup dominating_set
 * The class for the Schiermeyer exact dominating set algorithm.
 *
 */
class SchiermeyerDominatingSet : public ExactDominatingSet {
 public:
    using ExactDominatingSet::ExactDominatingSet;

    /**
     * Execute the Schiermeyer exact dominating set algorithm.
     */
    void run();

 private:
    bool find_small_MODS(const NetworKit::Graph &G);
    void find_big_MODS(const NetworKit::Graph &G);
    std::vector<bool> get_matching_MODS(const std::set<NetworKit::node> &S);
};

///////////

class ExhaustiveDominatingSet : public DominatingSet {
 public:
    using DominatingSet::DominatingSet;
    void run();
 private:
    std::vector<bool> recursiveDominatingSubset(std::vector<bool> &superset, int depth);
};

}  /* namespace Koala */
