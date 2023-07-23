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
        for (const auto &u : graph->nodeRange()) {
            std::set<NetworKit::node> neighborhood(
                graph->neighborRange(u).begin(), graph->neighborRange(u).end());
            neighborhood.insert(u);
            family.emplace_back(neighborhood);
        }
        std::vector<std::set<NetworKit::index>> occurences(family);
        auto set_cover_algorithm = SetCoverAlgorithm(family, occurences);
        set_cover_algorithm.run();
        auto set_cover = set_cover_algorithm.getSetCover();
        for (int i = 0; i < set_cover.size(); i++) {
            if (set_cover[i]) {
                dominating_set.insert(i);
            }
        }
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

    // TODO(kturowski): return solution or empty set
    bool find_small_MODS_recursive(
        const NetworKit::Graph &G, const std::vector<NetworKit::node> &V,
        NetworKit::index index, NetworKit::count size, std::set<NetworKit::node> &S);
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
    std::set<NetworKit::node> degree[2];

    std::set<NetworKit::node> find_big_MODS_recursive(NetworKit::Graph &G);
    std::set<NetworKit::node> find_MODS_for_minimum_degree_3(NetworKit::Graph &G);
    std::vector<NetworKit::node> move_to_solution(NetworKit::Graph &G, NetworKit::node vertex);
    void remove_from_solution(NetworKit::node vertex, const std::vector<NetworKit::node> &moved);
    bool forget_vertex(NetworKit::Graph &G, NetworKit::node vertex, bool is_required);
    void retrieve_vertex(
        NetworKit::Graph &G, NetworKit::node vertex, bool is_free, bool is_required);
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
    std::set<NetworKit::node> neighborhood;

    static NetworKit::Graph get_core_graph(
        const NetworKit::Graph &G, std::set<NetworKit::node> &free,
        std::set<NetworKit::node> &bound, std::set<NetworKit::node> &required);
    bool find_small_MODS(
        const NetworKit::Graph &G, const std::vector<NetworKit::node> &V);
    void find_big_MODS(
        const NetworKit::Graph &G, const std::vector<NetworKit::node> &V);
    std::vector<NetworKit::node> find_big_MODS_recursive(
        const NetworKit::Graph &G, const std::vector<NetworKit::node> &V, NetworKit::index index);
    std::vector<NetworKit::node> get_new_neighborhood(
        const NetworKit::Graph &G, NetworKit::node vertex);
    static std::vector<NetworKit::node> get_matching_MODS(
        const NetworKit::Graph &G, std::set<NetworKit::node> free,
        std::set<NetworKit::node> bound, std::set<NetworKit::node> required);
};

}  /* namespace Koala */
