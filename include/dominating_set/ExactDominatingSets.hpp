#pragma once

#include <set>

#include <dominating_set/MinimumDominatingSet.hpp>

template<typename BranchAndReduceMCS>
class BranchAndReduceMDS : public MinimumDominatingSet {
 public:
    using MinimumDominatingSet::MinimumDominatingSet;

    explicit BranchAndReduceMDS(const NetworKit::Graph &G) : MinimumDominatingSet(G) {}

    void run() override {
        std::vector<std::set<NetworKit::node>> family;
        G->forNodes([&family, this](NetworKit::node u) {
            std::set<NetworKit::node> neighborhood;
            neighborhood.insert(u);
            G->forNeighborsOf(u, [&neighborhood](NetworKit::node v) {
                neighborhood.insert(v);
            });
            family.emplace_back(neighborhood);
        });
        std::vector<std::set<NetworKit::index>> occurences(family);
        dominatingSet = BranchAndReduceMCS(family, occurences).run();
        hasRun = true;
    }
};

class FominKratschWoegingerMDS : public MinimumDominatingSet {
 public:
    using MinimumDominatingSet::MinimumDominatingSet;
    explicit FominKratschWoegingerMDS(const NetworKit::Graph &G);
    void run() override;
};

class ExhaustiveMDS : public MinimumDominatingSet {
 public:
    using MinimumDominatingSet::MinimumDominatingSet;
    explicit ExhaustiveMDS(const NetworKit::Graph &G);
    void run() override;
 private:
    std::vector<bool> recursiveDominatingSubset(std::vector<bool> &superset, int depth);
};

std::vector<NetworKit::node> joinFreeAndBounded(
    const std::set<NetworKit::node> &free,
    const std::set<NetworKit::node> &bounded);
bool isOptionalDominatingSet(
    const NetworKit::Graph &G,
    const std::vector<NetworKit::node> &choices,
    const std::set<NetworKit::node> &bounded);

class SchiermeyerMDS : public MinimumDominatingSet {
 public:
    using MinimumDominatingSet::MinimumDominatingSet;
    explicit SchiermeyerMDS(const NetworKit::Graph &G);
    void run() override;
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
