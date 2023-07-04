#include <cassert>

#include <set_cover/BranchAndReduceMSC.hpp>

bool isSubset(
        std::set<NetworKit::node> &subsetCandidate,
        std::set<NetworKit::node> &supersetCandidate) {
    if (subsetCandidate.size() > supersetCandidate.size()) {
        return false;
    }
    return std::all_of(
        subsetCandidate.begin(),
        subsetCandidate.end(),
        [&supersetCandidate](const auto& node) {
            return supersetCandidate.count(node) > 0;
        });
}

std::tuple<NetworKit::index, NetworKit::index> findSetInclusion(
        std::vector<std::set<NetworKit::count>> &sets) {
    for (NetworKit::count i = 0; i < sets.size(); i++) {
        auto &subsetCandidate = sets.at(i);
        if (subsetCandidate.empty()) {
            continue;
        }
        for (NetworKit::count j = 0; j < sets.size(); j++) {
            if (i == j) {
                continue;
            }
            auto &supersetCandidate = sets.at(j);
            if (isSubset(subsetCandidate, supersetCandidate)) {
                return {i, j};
            }
        }
    }
    return {NetworKit::none, NetworKit::none};
}

void excludeSet(
        NetworKit::index id,
        std::set<NetworKit::count> &excluded,
        std::vector<std::set<NetworKit::count>> &reversed) {
    for (auto &element : excluded) {
        reversed.at(element).erase(id);
    }
}

void includeSet(
        NetworKit::index id,
        std::set<NetworKit::count> &included,
        std::vector<std::set<NetworKit::count>> &reversed) {
    for (auto &element : included) {
        reversed.at(element).insert(id);
    }
}

RooijBodlaenderMSC::RooijBodlaenderMSC(
    std::vector<std::set<NetworKit::node>> &family,
    std::vector<std::set<NetworKit::index>> &occurences)
    : BranchAndReduceMSCImpl<true>(family, occurences) {}

NetworKit::index RooijBodlaenderMSC::findCountingRuleReductionSet() {
    for (NetworKit::index i = 0; i < family.size(); i++) {
        auto &candidate = family.at(i);
        int numberOfFrequency2elements = 0;
        std::set<NetworKit::node> notCoveredByCandidate;
        for (auto &element : candidate) {
            if (occurences.at(element).size() == 2) {
                numberOfFrequency2elements++;
                for (auto &occurence : occurences.at(element)) {
                    if (occurence == i) {
                        continue;
                    }
                    for (auto &e : family.at(occurence)) {
                        if (candidate.count(e) == 0) {
                            notCoveredByCandidate.insert(e);
                        }
                    }
                }
            }
        }
        if (notCoveredByCandidate.size() < numberOfFrequency2elements) {
            return i;
        }
    }
    return NetworKit::none;
}

NetworKit::index RooijBodlaenderMSC::find2Cardinality2FrequencySet() {
    for (NetworKit::index i = 0; i < family.size(); i++) {
        auto &candidate = family.at(i);
        if (candidate.size() != 2) {
            continue;
        }
        bool satifies = std::all_of(
            candidate.begin(),
            candidate.end(),
            [this] (const auto& element) {
                return occurences.at(element).size() == 2;
            });
        if (satifies) {
            return i;
        }
    }
    return NetworKit::none;
}

bool RooijBodlaenderMSC::reduce(std::vector<bool> & solution) {
    if (BranchAndReduceMSCImpl<true>::reduce(solution)) {
        return true;
    }
    // subsumption rule
    auto [subsetIndex2, supersetIndex2] = findSetInclusion(occurences);
    if (subsetIndex2 != NetworKit::none) {
        std::set<NetworKit::index> swapped;
        auto &superset = occurences.at(supersetIndex2);
        excludeSet(supersetIndex2, superset, family);
        superset.swap(swapped);
        std::vector<bool> subcover = run();
        superset.swap(swapped);
        includeSet(supersetIndex2, superset, family);
        solution = subcover;
        return true;
    }
    // counting rule
    NetworKit::index countingRuleIndex = findCountingRuleReductionSet();
    if (countingRuleIndex != NetworKit::none) {
        solution = forcedSetCover(countingRuleIndex);
        return true;
    }
    // size two set with frequency two elements rule
    NetworKit::index cardinalityFrequencyIndex = find2Cardinality2FrequencySet();
    if (cardinalityFrequencyIndex != NetworKit::none) {
        std::vector<NetworKit::index> indicesOfReplaced {cardinalityFrequencyIndex};
        for (auto &element : family.at(cardinalityFrequencyIndex)) {
            for (auto removed : occurences.at(element)) {
                if (removed == cardinalityFrequencyIndex) {
                    continue;
                }
                indicesOfReplaced.push_back(removed);
            }
        }
        std::set<NetworKit::node> replacement;
        auto& set1 = family.at(indicesOfReplaced[1]);
        auto& set2 = family.at(indicesOfReplaced[2]);
        std::set_union(
            set1.begin(),
            set1.end(),
            set2.begin(),
            set2.end(),
            std::inserter(replacement, replacement.begin()));
        for (auto &element : family.at(cardinalityFrequencyIndex)) {
            replacement.erase(element);
        }
        std::set<NetworKit::node> swapped[3]{ {}, replacement, {} };
        for (int j = 0; j < 3; j++) {
            excludeSet(indicesOfReplaced[j], family.at(indicesOfReplaced[j]), occurences);
        }
        includeSet(indicesOfReplaced[1], swapped[1], occurences);
        for (int j = 0; j < 3; j++) {
            swapped[j].swap(family.at(indicesOfReplaced[j]));
        }
        std::vector<bool> replacedCover = run();
        for (int j = 0; j < 3; j++) {
            swapped[j].swap(family.at(indicesOfReplaced[j]));
        }
        excludeSet(indicesOfReplaced[1], swapped[1], occurences);
        for (int j = 0; j < 3; j++) {
            includeSet(indicesOfReplaced[j], family.at(indicesOfReplaced[j]), occurences);
        }
        if (replacedCover.at(indicesOfReplaced[1])) {
            replacedCover.at(indicesOfReplaced[2]) = true;
        } else {
            replacedCover.at(indicesOfReplaced[0]) = true;
        }
        solution = replacedCover;
        return true;
    }
    return false;
}
