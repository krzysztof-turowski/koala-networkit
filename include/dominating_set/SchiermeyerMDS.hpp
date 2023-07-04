#pragma once

#include <set>

#include <dominating_set/MinimumDominatingSet.hpp>

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
