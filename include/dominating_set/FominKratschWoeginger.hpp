#pragma once

#include <set>

#include <dominating_set/MinimumDominatingSet.hpp>

class FominKratschWoegingerMDS : public MinimumDominatingSet {
 public:
    using MinimumDominatingSet::MinimumDominatingSet;
    explicit FominKratschWoegingerMDS(const NetworKit::Graph &G);
    void run() override;
};
