#pragma once

#include <set>

#include <dominating_set/MinimumDominatingSet.hpp>

class ExhaustiveMDS : public MinimumDominatingSet {
 public:
    using MinimumDominatingSet::MinimumDominatingSet;
    explicit ExhaustiveMDS(const NetworKit::Graph &G);
    void run() override;
 private:
    std::vector<bool> recursiveDominatingSubset(std::vector<bool> &superset, int depth);
};
