#include <dominating_set/MinimumDominatingSet.hpp>

namespace Koala {

MinimumDominatingSet::MinimumDominatingSet(const NetworKit::Graph &G) : G(&G) {}

bool MinimumDominatingSet::isDominating(const std::vector<bool> &dominating_set) {
    std::vector<bool> dominated(dominating_set.size());
    G->forNodes([&dominated, &dominating_set, this](NetworKit::node u) {
        if (dominating_set[u]) {
            dominated[u] = true;
            G->forNeighborsOf(u, [&dominated, u](NetworKit::node v) {
                dominated[v] = true;
            });
        }
    });
    return std::all_of(
        dominated.begin(),
        dominated.end(),
        [](bool is_dominated) {
            return is_dominated;
        });
}

int MinimumDominatingSet::dominatingSetSize(const std::vector<bool> &set) {
    return std::count(set.begin(), set.end(), true);
}

std::vector<bool> &smallerCardinalitySet(std::vector<bool> &lhs, std::vector<bool> &rhs) {
    if (std::count(lhs.begin(), lhs.end(), true) < std::count(rhs.begin(), rhs.end(), true)) {
        return lhs;
    } else {
        return rhs;
    }
}
}  /* namespace Koala */
