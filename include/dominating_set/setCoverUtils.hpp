#pragma once

#include <tuple>
#include <set>

#include <dominating_set/MinimumDominatingSet.hpp>

std::tuple<NetworKit::index, NetworKit::index> findSetInculsion(
    std::vector<std::set<NetworKit::count>> &sets);
void excludeSet(
    NetworKit::index id, std::set<NetworKit::node> &excluded,
    std::vector<std::set<NetworKit::index>> &reversed);
void includeSet(
    NetworKit::index id, std::set<NetworKit::node> &included,
    std::vector<std::set<NetworKit::index>> &reversed);

bool isSubset(
        std::set<NetworKit::node> &subsetCandidate,
        std::set<NetworKit::node> &supersetCandidate) {
    if (subsetCandidate.size() > supersetCandidate.size()) {
        return false;
    }
    for (auto &node : subsetCandidate) {
        if (supersetCandidate.count(node) == 0) {
            return false;
        }
    }
    return true;
}
