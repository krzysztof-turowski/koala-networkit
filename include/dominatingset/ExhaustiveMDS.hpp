#ifndef EXHAUSTIVE_MDS_HPP_
#define EXHAUSTIVE_MDS_HPP_

#include <set>
#include <dominatingset/MinimumDominatingSet.hpp>

class ExhaustiveMDS : public MinimumDominatingSet {
public:
    using MinimumDominatingSet::MinimumDominatingSet;
    ExhaustiveMDS(const NetworKit::Graph &G);
    void run() override;
private:
    std::vector<bool> reccursiveDominatingSubset(std::vector<bool> &superset, int depth);
};

#endif /* EXHAUSTIVE_MDS_HPP_ */
