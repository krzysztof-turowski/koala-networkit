#pragma once

#include <set>

#include <dominating_set/MinimumDominatingSet.hpp>

class RooijBodlaenderMDS : public MinimumDominatingSet {
 public:
    using MinimumDominatingSet::MinimumDominatingSet;
    explicit RooijBodlaenderMDS(const NetworKit::Graph &G);
    void run() override;
};
