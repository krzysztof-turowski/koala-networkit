#include <dominating_set/ExactDominatingSet.hpp>

namespace Koala {

std::vector<bool> &smallerCardinalitySet(std::vector<bool> &lhs, std::vector<bool> &rhs) {
    if (std::count(lhs.begin(), lhs.end(), true) < std::count(rhs.begin(), rhs.end(), true)) {
        return lhs;
    } else {
        return rhs;
    }
}

void ExhaustiveMDS::run() {
    std::vector<bool> all_vertices;
    graph->forNodes([&all_vertices, this](NetworKit::node u) {
        all_vertices.emplace_back(true);
    });

    dominating_set = recursiveDominatingSubset(all_vertices, 0);
    hasRun = true;
}

std::vector<bool> ExhaustiveMDS::recursiveDominatingSubset(std::vector<bool> &superset, int depth) {
    if (depth == superset.size()) {
        return superset;
    }
    superset[depth] = false;

    std::vector<bool> dominated(superset.size());
    graph->forNodes([&dominated, &superset, this](NetworKit::node u) {
        if (superset[u]) {
            dominated[u] = true;
            graph->forNeighborsOf(u, [&dominated](NetworKit::node v) {
                dominated[v] = true;
            });
        }
    });
    bool isDominating = std::all_of(dominated.begin(), dominated.end(), [](bool d) { return d; });
    if (isDominating) {
        std::vector<bool> excluded = recursiveDominatingSubset(superset, depth + 1);
        superset[depth] = true;
        std::vector<bool> included = recursiveDominatingSubset(superset, depth + 1);
        return smallerCardinalitySet(excluded, included);
    } else {
        superset[depth] = true;
        return recursiveDominatingSubset(superset, depth + 1);
    }
}
}  /* namespace Koala */
