#include <cassert>

#include <dominating_set/BranchAndReduceMSC.hpp>
#include <dominating_set/RooijBodlaenderMDS.hpp>

NetworKit::index findCountingRuleReductionSet(
    std::vector<std::set<NetworKit::node>> &family,
    std::vector<std::set<NetworKit::index>> &occurences);
NetworKit::index find2Cardinality2FrequencySet(
    std::vector<std::set<NetworKit::node>> &family,
    std::vector<std::set<NetworKit::index>> &occurences);

RooijBodlaenderMSC::RooijBodlaenderMSC(
    std::vector<std::set<NetworKit::node>> &family,
    std::vector<std::set<NetworKit::index>> &occurences)
    : BranchAndReduceMSCImpl<true>(family, occurences) {}

bool RooijBodlaenderMSC::reduce(std::vector<bool> & solution) {
    MinimumDominatingSet::specialCounter2++;
    if (BranchAndReduceMSCImpl<true>::reduce(solution)) {
        return true;
    }
    // subsumption rule
    auto [subsetIndex2, supersetIndex2] = findSetInculsion(occurences);
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
    NetworKit::index countingRuleIndex = findCountingRuleReductionSet(family, occurences);
    if (countingRuleIndex != NetworKit::none) {
        solution = forcedSetCover(countingRuleIndex);
        return true;
    }
    // size two set with frequency two elements rule
    NetworKit::index cardinalityFrequencyIndex = find2Cardinality2FrequencySet(family, occurences);
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
