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

template<typename SetCoverAlgorithm>
class BranchAndReduceDominatingSet : public DominatingSet {
 public:
    using DominatingSet::DominatingSet;

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

class FominKratschWoegingerMDS : public DominatingSet {
 public:
    using DominatingSet::DominatingSet;
    void run();
};

class ExhaustiveMDS : public DominatingSet {
 public:
    using DominatingSet::DominatingSet;
    void run();
 private:
    std::vector<bool> recursiveDominatingSubset(std::vector<bool> &superset, int depth);
};

class SchiermeyerMDS : public DominatingSet {
 public:
    using DominatingSet::DominatingSet;
    void run();
};

class SizedChoiceSearcher {
    const std::function<bool(const std::vector<NetworKit::node>)> &verifier;
    const std::vector<NetworKit::node> &possibilities;
    int size;
    bool recursive(std::vector<NetworKit::node> &choices, NetworKit::node decideOn, int left);
 public:
    SizedChoiceSearcher(
        const std::function<bool(const std::vector<NetworKit::node>)> &verifier,
        const std::vector<NetworKit::node> &possibilities,
        int size)
        : verifier(verifier), possibilities(possibilities), size(size) {}
    std::tuple<bool, std::vector<NetworKit::node>> search();
};

}  /* namespace Koala */
