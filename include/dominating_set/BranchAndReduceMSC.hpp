#pragma once

#include <set>
#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/max_cardinality_matching.hpp>

#include <networkit/graph/Graph.hpp>

#include <dominating_set/MinimumDominatingSet.hpp>

template<bool useEdgeCover>
class BranchAndReduceMSCImpl {
 protected:
    std::vector<std::set<NetworKit::node>> &family;
    std::vector<std::set<NetworKit::index>> &occurences;
    virtual bool reduce(std::vector<bool> & solution) {
        // unique element rule
        NetworKit::index forcedIndex = findUniqueOccurenceSet();
        if (forcedIndex != NetworKit::none) {
            solution = forcedSetCover(forcedIndex);
            return true;
        }
        // subset rule
        auto [subsetIndex1, supersetIndex1] = findSetInculsion(family);
        if (subsetIndex1 != NetworKit::none) {
            solution = discardedSetCover(subsetIndex1);
            return true;
        }
        return false;
    }

    std::vector<bool> discardedSetCover(NetworKit::index discardedIndex) {
        std::set<NetworKit::node> swapped;
        excludeSet(discardedIndex, family.at(discardedIndex), occurences);
        family.at(discardedIndex).swap(swapped);
        std::vector<bool> largestExcluded = run();
        family.at(discardedIndex).swap(swapped);
        includeSet(discardedIndex, family.at(discardedIndex), occurences);
        return largestExcluded;
    }

    std::vector<bool> forcedSetCover(NetworKit::index forcedIndex) {
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
        std::vector<bool> subcover = run();
        for (auto &element : removed) {
            family.at(element.first).insert(element.second);
            occurences.at(element.second).insert(element.first);
        }
        subcover.at(forcedIndex) = true;
        return subcover;
    }

    std::vector<bool> blossom() {
        typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> boost_graph_t;
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
    NetworKit::index findUniqueOccurenceSet() {
        for (auto &elementOccurences : occurences) {
            if (elementOccurences.size() == 1) {
                NetworKit::node forcedIndex = *(elementOccurences.begin());
                return forcedIndex;
            }
        }
        return NetworKit::none;
    }

    bool reducesToMatching() {
        for (auto &familyElement : family) {
            if (familyElement.size() > 2) {
                return false;
            }
        }
        return true;
    }

 public:
    BranchAndReduceMSCImpl(
        std::vector<std::set<NetworKit::node>> &family,
        std::vector<std::set<NetworKit::index>> &occurences)
        : family(family), occurences(occurences) {}
    std::vector<bool> run() {
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
        std::vector<bool> solutionFromReduce;
        if (reduce(solutionFromReduce)) {
            return solutionFromReduce;
        }
        if (useEdgeCover && reducesToMatching()) {
            return blossom();
        }
        // branching
        auto iteratorOfLargest = std::max_element(
            family.begin(),
            family.end(),
            [](const auto& lhs, const auto& rhs) {
                return lhs.size() < rhs.size();
            });
        NetworKit::index indexOfLargest = iteratorOfLargest - family.begin();
        std::vector<bool> largestIncluded = forcedSetCover(indexOfLargest);
        std::vector<bool> largestExcluded = discardedSetCover(indexOfLargest);

        return smallerCardinalitySet(largestIncluded, largestExcluded);
    }
};

typedef BranchAndReduceMSCImpl<false> GrandoniMSC;
typedef BranchAndReduceMSCImpl<true> FominGrandoniKratschMSC;

class RooijBodlaenderMSC : public BranchAndReduceMSCImpl<true> {
    NetworKit::index findCountingRuleReductionSet();
    NetworKit::index find2Cardinality2FrequencySet();
 protected:
    bool reduce(std::vector<bool> &solution);
 public:
    RooijBodlaenderMSC(
        std::vector<std::set<NetworKit::node>> &family,
        std::vector<std::set<NetworKit::index>> &occurences);
};

std::tuple<NetworKit::index, NetworKit::index> findSetInculsion(
    std::vector<std::set<NetworKit::count>> &sets);
void excludeSet(
    NetworKit::index id, std::set<NetworKit::node> &excluded,
    std::vector<std::set<NetworKit::index>> &reversed);
void includeSet(
    NetworKit::index id, std::set<NetworKit::node> &included,
    std::vector<std::set<NetworKit::index>> &reversed);