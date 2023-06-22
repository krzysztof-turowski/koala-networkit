#ifndef SCHIERMEYER_MDS_HPP_
#define SCHIERMEYER_MDS_HPP_

#include <set>
#include <dominatingset/MinimumDominatingSet.hpp>

std::vector<NetworKit::node> joinFreeAndBounded(const std::set<NetworKit::node> &free, const std::set<NetworKit::node> &bounded);
bool isOptionalDominatingSet(const NetworKit::Graph &G, const std::vector<NetworKit::node> &choices, const std::set<NetworKit::node> &bounded);
bool recursiveSizedChoiceSearch(const std::function<bool(const std::vector<NetworKit::node>)>&  verifier, const std::vector<NetworKit::node> &possibilities, std::vector<NetworKit::node> &choices, const NetworKit::node decideOn, const int left);

class SchiermeyerMDS : public MinimumDominatingSet {
public:
    using MinimumDominatingSet::MinimumDominatingSet;
    SchiermeyerMDS(const NetworKit::Graph &G);
    void run() override;
};

#endif /* SCHIERMEYER_MDS_HPP_ */
