#ifndef BRANCH_AND_REDUCE_MDS_HPP_
#define BRANCH_AND_REDUCE_MDS_HPP_

#include <set>
#include <dominatingset/MinimumDominatingSet.hpp>

template<typename BranchAndReduceMCS>
class BranchAndReduceMDS : public MinimumDominatingSet {
public:
    using MinimumDominatingSet::MinimumDominatingSet;

    BranchAndReduceMDS(const NetworKit::Graph &G) : MinimumDominatingSet(G) {}

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

#endif /* BRANCH_AND_REDUCE_MDS_HPP_ */
