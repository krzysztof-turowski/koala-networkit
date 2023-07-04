#include <dominating_set/ExhaustiveMDS.hpp>

ExhaustiveMDS::ExhaustiveMDS(const NetworKit::Graph &G) : MinimumDominatingSet(G) {}

void ExhaustiveMDS::run() {
    std::vector<bool> all_vertices;
    G->forNodes([&all_vertices, this](NetworKit::node u) {
        all_vertices.emplace_back(true);
    });

    dominatingSet = recursiveDominatingSubset(all_vertices, 0);

    hasRun = true;
}

std::vector<bool> ExhaustiveMDS::recursiveDominatingSubset(std::vector<bool> &superset, int depth) {
    if (depth == superset.size()) {
        return superset;
    }
    superset[depth] = false;
    if (isDominating(superset)) {
        std::vector<bool> excluded = recursiveDominatingSubset(superset, depth + 1);
        superset[depth] = true;
        std::vector<bool> included = recursiveDominatingSubset(superset, depth + 1);
        return smallerCardinalitySet(excluded, included);
    } else {
        superset[depth] = true;
        return recursiveDominatingSubset(superset, depth + 1);
    }
}
