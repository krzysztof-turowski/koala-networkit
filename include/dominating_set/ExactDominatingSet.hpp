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
        set_cover_algorithm.getSetCover().swap(dominating_set);
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
    std::set<NetworKit::node> free, bound, required;
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
    std::set<NetworKit::node> degree_one, degree_two;
    std::vector<std::set<NetworKit::node>> neighborhood;
    std::vector<NetworKit::node> possibilities;

    std::vector<bool> find_big_MODS_recursive();
    std::vector<bool> find_MODS_for_minimum_degree_3();
    std::tuple<std::vector<NetworKit::node>, bool> add_to_solution(NetworKit::node vertex);
    void remove_from_solution(
        NetworKit::node vertex, std::vector<NetworKit::node> &moved, bool is_free);
    bool forget_vertex(NetworKit::node vertex);
    void retrieve_vertex(NetworKit::node vertex, bool is_free);
    void on_degree_decrement(NetworKit::node vertex);
    void on_degree_increment(NetworKit::node vertex);
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

    std::set<NetworKit::node> choices, neighborhood;
    std::vector<NetworKit::node> possibilities;

    void find_big_MODS_recursive(const NetworKit::Graph &G, NetworKit::count index);
    std::set<NetworKit::node> get_new_neighborhood(
        const NetworKit::Graph &G, NetworKit::node vertex);

    std::vector<bool> get_matching_MODS(
        const NetworKit::Graph &G, const std::set<NetworKit::node> &required);
    NetworKit::Graph core(
        const NetworKit::Graph &G, std::set<NetworKit::node> &free,
        std::set<NetworKit::node> &bound, std::set<NetworKit::node> &required);
};

}  /* namespace Koala */
