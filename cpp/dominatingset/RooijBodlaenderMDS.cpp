#include <dominatingset/RooijBodlaenderMDS.hpp>
#include <cassert>
#include <boost/graph/max_cardinality_matching.hpp>
#include <boost/graph/adjacency_list.hpp>

bool isSubset(std::set<NetworKit::node> &subsetCandidate, std::set<NetworKit::node> &supersetCandidate);
std::vector<bool> minimumSetCover(std::vector<std::set<NetworKit::node>> &family, std::vector<std::set<NetworKit::index>> &occurences);
NetworKit::index findUniqueOccurenceSet(std::vector<std::set<NetworKit::index>> &occurences);
NetworKit::index findCountingRuleReductionSet(std::vector<std::set<NetworKit::node>> &family, std::vector<std::set<NetworKit::index>> &occurences);
NetworKit::index find2Cardinality2FrequencySet(std::vector<std::set<NetworKit::node>> &family, std::vector<std::set<NetworKit::index>> &occurences);
std::tuple<NetworKit::index, NetworKit::index> findSetInculsion(std::vector<std::set<NetworKit::count>> &sets);
bool reducesToMatching(std::vector<std::set<NetworKit::node>> &family);
std::vector<bool> discardedSetCover(NetworKit::index discardedIndex, std::vector<std::set<NetworKit::node>> &family, std::vector<std::set<NetworKit::index>> &occurences);
std::vector<bool> forcedSetCover(NetworKit::index forcedIndex, std::vector<std::set<NetworKit::node>> &family, std::vector<std::set<NetworKit::index>> &occurences);
void excludeSet(NetworKit::index id, std::set<NetworKit::node> &excluded, std::vector<std::set<NetworKit::index>> &occurences);
void includeSet(NetworKit::index id, std::set<NetworKit::node> &included, std::vector<std::set<NetworKit::index>> &occurences);
std::vector<bool> blossom(std::vector<std::set<NetworKit::node>> &family, std::vector<std::set<NetworKit::index>> &occurences);

template <typename T>
std::set<T> setUnion(const std::set<T>& first, const std::set<T>& second)
{
    std::set<T> result = first;
    result.insert(second.begin(), second.end());
    return result;
}

RooijBodlaenderMDS::RooijBodlaenderMDS(const NetworKit::Graph &G) : MinimumDominatingSet(G) {}

void RooijBodlaenderMDS::run() {
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
    dominatingSet = minimumSetCover(family, occurences);
    hasRun = true;
}

std::vector<bool> minimumSetCover(std::vector<std::set<NetworKit::node>> &family, std::vector<std::set<NetworKit::index>> &occurences) {
    bool emptyInstance = true;
    for (auto &familyElement : family) {
        if (familyElement.size() != 0) {
            emptyInstance = false;
            break;
        }
    }
    if (emptyInstance) {
        return std::vector<bool>(family.size());
    }
    MinimumDominatingSet::specialCounter1++;
    
    // unique element rule
    NetworKit::index forcedIndex = findUniqueOccurenceSet(occurences);
    if (forcedIndex != NetworKit::none) {
        return forcedSetCover(forcedIndex, family, occurences);
    }
    
    // subset rule
    auto [subsetIndex1, supersetIndex1] = findSetInculsion(family);
    if (subsetIndex1 != NetworKit::none) {
        return discardedSetCover(subsetIndex1, family, occurences);
    }

    // subsumption rule
    auto [subsetIndex2, supersetIndex2] = findSetInculsion(occurences);
    if (subsetIndex2 != NetworKit::none) {
        std::set<NetworKit::index> swapped;
        auto &superset = occurences.at(supersetIndex2);
        excludeSet(supersetIndex2, superset, family);
        superset.swap(swapped);
        std::vector<bool> subcover = minimumSetCover(family, occurences);
        superset.swap(swapped);
        includeSet(supersetIndex2, superset, family);
        return subcover;
    }

    // counting rule
    NetworKit::index countingRuleIndex = findCountingRuleReductionSet(family, occurences);
    if (countingRuleIndex != NetworKit::none) {
        return forcedSetCover(countingRuleIndex, family, occurences);
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
        std::set<NetworKit::node> replacement = setUnion<NetworKit::node>(family.at(indicesOfReplaced[1]), family.at(indicesOfReplaced[2]));
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
        std::vector<bool> replacedCover = minimumSetCover(family, occurences);
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
        return replacedCover;
    }

    if (reducesToMatching(family)) {
        return blossom(family, occurences);
    }

    //branching
    auto iteratorOfLargest = std::max_element(family.begin(), family.end(), [](const std::set<NetworKit::node>& lhs, const std::set<NetworKit::node>& rhs) {
        return lhs.size() < rhs.size();
    });
    NetworKit::index indexOfLargest = iteratorOfLargest - family.begin();
    std::vector<bool> largestIncluded = forcedSetCover(indexOfLargest, family, occurences);
    std::vector<bool> largestExcluded = discardedSetCover(indexOfLargest, family, occurences);

    return smallerCardinalitySet(largestIncluded, largestExcluded);
}

NetworKit::index findUniqueOccurenceSet(std::vector<std::set<NetworKit::index>> &occurences) {
    for (auto &elementOccurences : occurences) {
        if (elementOccurences.size() == 1) {
            NetworKit::node forcedIndex = *(elementOccurences.begin());
            return forcedIndex;
        }
    }
    return NetworKit::none;
}

NetworKit::index findCountingRuleReductionSet(std::vector<std::set<NetworKit::node>> &family, std::vector<std::set<NetworKit::index>> &occurences) {
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

NetworKit::index find2Cardinality2FrequencySet(std::vector<std::set<NetworKit::node>> &family, std::vector<std::set<NetworKit::index>> &occurences) {
    for (NetworKit::index i = 0; i < family.size(); i++) {
        auto &candidate = family.at(i);
        if (candidate.size() != 2) {
            continue;
        }
        bool satifies = true;
        for (auto &element : candidate) {
            if (occurences.at(element).size() != 2) {
                satifies = false;
                break;
            }
        }
        if (satifies) {
            return i;
        }
    }
    return NetworKit::none;
}

std::tuple<NetworKit::index, NetworKit::index> findSetInculsion(std::vector<std::set<NetworKit::count>> &sets) {
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

bool reducesToMatching(std::vector<std::set<NetworKit::node>> &family) {
    for (auto &familyElement : family) {
        if (familyElement.size() > 2) {
            return false;
        }
    }
    return true;
}

std::vector<bool> discardedSetCover(NetworKit::index discardedIndex, std::vector<std::set<NetworKit::node>> &family, std::vector<std::set<NetworKit::index>> &occurences) {
    std::set<NetworKit::node> swapped;
    excludeSet(discardedIndex, family.at(discardedIndex), occurences);
    family.at(discardedIndex).swap(swapped);
    std::vector<bool> largestExcluded = minimumSetCover(family, occurences);
    family.at(discardedIndex).swap(swapped);
    includeSet(discardedIndex, family.at(discardedIndex), occurences);
    return largestExcluded;
}

std::vector<bool> forcedSetCover(NetworKit::index forcedIndex, std::vector<std::set<NetworKit::node>> &family, std::vector<std::set<NetworKit::index>> &occurences) {
    std::set<NetworKit::node> &forced = family.at(forcedIndex);
    std::vector<std::pair<NetworKit::index, NetworKit::node>> removed;
    for (auto &element : forced) {
        for (auto &elementOccurence : occurences.at(element)) {
            removed.emplace_back(elementOccurence, element);
        }
        occurences.at(element) = std::set<NetworKit::index>();
    }
    for (auto &element : removed) {
        family.at(element.first).erase(element.second);
    }
    std::vector<bool> subcover = minimumSetCover(family, occurences);
    for (auto &element : removed) {
        family.at(element.first).insert(element.second);
        occurences.at(element.second).insert(element.first);
    }
    subcover.at(forcedIndex) = true;
    return subcover;
}

bool isSubset(std::set<NetworKit::node> &subsetCandidate, std::set<NetworKit::node> &supersetCandidate) {
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

void excludeSet(NetworKit::index id, std::set<NetworKit::count> &excluded, std::vector<std::set<NetworKit::count>> &reversed) {
    for (auto &element : excluded) {
        reversed.at(element).erase(id);
    }
}

void includeSet(NetworKit::index id, std::set<NetworKit::count> &included, std::vector<std::set<NetworKit::count>> &reversed) {
    for (auto &element : included) {
        reversed.at(element).insert(id);
    }
}

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> boost_graph_t;
std::vector<bool> blossom(std::vector<std::set<NetworKit::node>> &family, std::vector<std::set<NetworKit::index>> &occurences) {
    std::vector<bool> cover(family.size());
    std::vector<bool> dominated(occurences.size());
    boost_graph_t boost_graph(occurences.size());
    std::vector<boost::graph_traits<boost_graph_t>::vertex_descriptor> mate(occurences.size());
    for (auto &element : family) {
        if (element.size() == 0) {
            continue;
        }
        std::set<NetworKit::node>::iterator it = element.begin();
        NetworKit::node v1 = *(it);
        NetworKit::node v2 = *(++it);
        boost::add_edge(v1, v2, boost_graph);
    }

    bool success = boost::checked_edmonds_maximum_cardinality_matching(boost_graph, &mate[0]);
    assert(success);

    for (NetworKit::index i = 0; i < family.size(); i++) {
        auto &element = family.at(i);
        if (element.size() == 0) {
            continue;
        }
        std::set<NetworKit::node>::iterator it = element.begin();
        NetworKit::node v1 = *(it);
        NetworKit::node v2 = *(++it);
        if (mate[v1] == v2) {
            cover[i] = true;
            dominated[v1] = true;
            dominated[v2] = true;
        }
    }
    for (NetworKit::node i = 0; i < occurences.size(); i++) {
        if (!dominated[i] && occurences.at(i).size() != 0) {
            cover[*occurences.at(i).begin()] = true;
        }
    }
    return cover;
}
