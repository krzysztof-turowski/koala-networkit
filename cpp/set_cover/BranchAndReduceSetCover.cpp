/*
 * BranchAndReduceSetCover.cpp
 *
 *  Created on: 01.07.2023
 *      Author: Piotr Kubaty
 */

#include <cassert>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/max_cardinality_matching.hpp>

#include <set_cover/BranchAndReduceSetCover.hpp>

namespace Koala {

bool is_subset(
        std::set<NetworKit::node> &subset, std::set<NetworKit::node> &superset) {
    if (subset.size() > superset.size()) {
        return false;
    }
    return std::all_of(
        subset.begin(), subset.end(), [&](const auto& v) { return superset.count(v) > 0; });
}

std::tuple<NetworKit::index, NetworKit::index> find_set_inclusion(
        std::vector<std::set<NetworKit::count>> &S) {
    for (NetworKit::count i = 0; i < S.size(); i++) {
        auto &subset = S.at(i);
        if (subset.empty()) {
            continue;
        }
        for (NetworKit::count j = 0; j < S.size(); j++) {
            if (i == j) {
                continue;
            }
            auto &superset = S.at(j);
            if (is_subset(subset, superset)) {
                return {i, j};
            }
        }
    }
    return {NetworKit::none, NetworKit::none};
}

void exclude_set(
        NetworKit::index id,
        std::set<NetworKit::count> &excluded,
        std::vector<std::set<NetworKit::count>> &reversed) {
    for (auto &element : excluded) {
        reversed.at(element).erase(id);
    }
}

void include_set(
        NetworKit::index id,
        std::set<NetworKit::count> &included,
        std::vector<std::set<NetworKit::count>> &reversed) {
    for (auto &element : included) {
        reversed.at(element).insert(id);
    }
}

BranchAndReduceSetCover::BranchAndReduceSetCover(
        std::vector<std::set<NetworKit::node>> &family,
        std::vector<std::set<NetworKit::index>> &occurences)
    : family(family), occurences(occurences) { }

std::vector<bool> BranchAndReduceSetCover::getSetCover() const {
    assureFinished();
    return set_cover;
}

void BranchAndReduceSetCover::check() const {
    assureFinished();
    // TODO
}

void BranchAndReduceSetCover::run() {
    hasRun = true;
    recurse().swap(set_cover);
}

std::vector<bool> BranchAndReduceSetCover::recurse() {
    if (std::all_of(family.begin(), family.end(), [](const auto& e) { return e.empty(); })) {
        return std::vector<bool>(family.size());
    }
    if (reduce()) {
        return set_cover;
    }
    if (reduce_matching()) {
        return set_cover;
    }
    // branching
    NetworKit::index largest = std::max_element(
        family.begin(), family.end(),
        [](const auto& a, const auto& b) { return a.size() < b.size(); }) - family.begin();
    auto included(forced_set_cover(largest)), excluded(discarded_set_cover(largest));
    int included_size = std::count(included.begin(), included.end(), true);
    int excluded_size = std::count(excluded.begin(), excluded.end(), true);
    set_cover.swap(included_size < excluded_size ? included : excluded);
    return set_cover;
}

bool BranchAndReduceSetCover::reduce() {
    // unique element rule
    NetworKit::index forced_index = find_unique_occurence_set();
    if (forced_index != NetworKit::none) {
        forced_set_cover(forced_index).swap(set_cover);
        return true;
    }
    // subset rule
    auto [subset_index, superset_index] = find_set_inclusion(family);
    if (subset_index != NetworKit::none) {
        discarded_set_cover(subset_index).swap(set_cover);
        return true;
    }
    return false;
}

NetworKit::index BranchAndReduceSetCover::find_unique_occurence_set() {
    for (auto &occurence : occurences) {
        if (occurence.size() == 1) {
            return *(occurence.begin());
        }
    }
    return NetworKit::none;
}

std::vector<bool> BranchAndReduceSetCover::forced_set_cover(NetworKit::index index) {
    auto &forced = family.at(index);
    std::vector<std::pair<NetworKit::index, NetworKit::node>> removed;
    for (auto &element : forced) {
        for (auto &occurence : occurences.at(element)) {
            removed.emplace_back(occurence, element);
        }
        occurences.at(element).clear();
    }
    for (auto &element : removed) {
        family.at(element.first).erase(element.second);
    }
    std::vector<bool> solution;
    recurse().swap(solution);
    for (auto &element : removed) {
        family.at(element.first).insert(element.second);
        occurences.at(element.second).insert(element.first);
    }
    solution.at(index) = true;
    return solution;
}

std::vector<bool> BranchAndReduceSetCover::discarded_set_cover(NetworKit::index index) {
    std::set<NetworKit::node> swapped;
    exclude_set(index, family.at(index), occurences);
    family.at(index).swap(swapped);
    std::vector<bool> solution;
    recurse().swap(solution);
    family.at(index).swap(swapped);
    include_set(index, family.at(index), occurences);
    return solution;
}

bool GrandoniSetCover::reduce_matching() {
    return false;
}

bool FominGrandoniKratschSetCover::reduce_matching() {
    if (std::any_of(family.begin(), family.end(), [](const auto& e) { return e.size() > 2; })) {
        return false;
    }

    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> boost_graph_t;
    boost_graph_t boost_graph(occurences.size());
    std::vector<boost::graph_traits<boost_graph_t>::vertex_descriptor> mate(occurences.size());
    for (auto &element : family) {
        if (element.empty()) {
            continue;
        }
        std::set<NetworKit::node>::iterator it = element.begin();
        boost::add_edge(*it, *std::next(it), boost_graph);
    }
    assert(boost::checked_edmonds_maximum_cardinality_matching(boost_graph, &mate[0]));

    set_cover.resize(family.size());
    std::vector<bool> dominated(occurences.size());
    for (NetworKit::index i = 0; i < family.size(); i++) {
        auto &element = family.at(i);
        if (element.empty()) {
            continue;
        }
        std::set<NetworKit::node>::iterator it = element.begin();
        NetworKit::node v1 = *it, v2 = *std::next(it);
        if (mate[v1] == v2) {
            set_cover[i] = dominated[v1] = dominated[v2] = true;
        }
    }
    for (NetworKit::node i = 0; i < occurences.size(); i++) {
        if (!dominated[i] && !occurences.at(i).empty()) {
            set_cover[*occurences.at(i).begin()] = true;
        }
    }
    return true;
}

bool RooijBodlaenderSetCover::reduce() {
    if (FominGrandoniKratschSetCover::reduce()) {
        return true;
    }
    // subsumption rule
    auto [subset_index, superset_index] = find_set_inclusion(occurences);
    if (subset_index != NetworKit::none) {
        std::set<NetworKit::index> temp;
        auto &superset = occurences.at(superset_index);
        exclude_set(superset_index, superset, family);
        superset.swap(temp);
        recurse().swap(set_cover);
        superset.swap(temp);
        include_set(superset_index, superset, family);
        return true;
    }
    // counting rule
    NetworKit::index counting_rule_index = find_counting_rule_reduction_set();
    if (counting_rule_index != NetworKit::none) {
        forced_set_cover(counting_rule_index).swap(set_cover);
        return true;
    }
    // size two set with frequency two elements rule
    NetworKit::index cardinality_frequency_index = find_cardinality_frequency_set();
    if (cardinality_frequency_index == NetworKit::none) {
        return false;
    }
    std::vector<NetworKit::index> indices{cardinality_frequency_index};
    for (auto &element : family.at(cardinality_frequency_index)) {
        for (auto removed : occurences.at(element)) {
            if (removed == cardinality_frequency_index) {
                continue;
            }
            indices.push_back(removed);
        }
    }
    std::set<NetworKit::node> replacement;
    auto &A = family.at(indices[1]), &B = family.at(indices[2]);
    std::set_union(
        A.begin(), A.end(), B.begin(), B.end(),
        std::inserter(replacement, replacement.begin()));
    for (auto &element : family.at(cardinality_frequency_index)) {
        replacement.erase(element);
    }
    std::set<NetworKit::node> swapped[3]{ {}, replacement, {} };
    for (int j = 0; j < 3; j++) {
        exclude_set(indices[j], family.at(indices[j]), occurences);
    }
    include_set(indices[1], swapped[1], occurences);
    for (int j = 0; j < 3; j++) {
        swapped[j].swap(family.at(indices[j]));
    }
    recurse().swap(set_cover);
    for (int j = 0; j < 3; j++) {
        swapped[j].swap(family.at(indices[j]));
    }
    exclude_set(indices[1], swapped[1], occurences);
    for (int j = 0; j < 3; j++) {
        include_set(indices[j], family.at(indices[j]), occurences);
    }
    if (set_cover.at(indices[1])) {
        set_cover.at(indices[2]) = true;
    } else {
        set_cover.at(indices[0]) = true;
    }
    return true;
}

NetworKit::index RooijBodlaenderSetCover::find_counting_rule_reduction_set() {
    for (NetworKit::index i = 0; i < family.size(); i++) {
        auto &candidate = family.at(i);
        int count_frequency_two = 0;
        std::set<NetworKit::node> not_covered;
        for (auto &element : candidate) {
            if (occurences.at(element).size() == 2) {
                count_frequency_two++;
                for (auto &occurence : occurences.at(element)) {
                    if (occurence == i) {
                        continue;
                    }
                    for (auto &e : family.at(occurence)) {
                        if (!candidate.count(e)) {
                            not_covered.insert(e);
                        }
                    }
                }
            }
        }
        if (not_covered.size() < count_frequency_two) {
            return i;
        }
    }
    return NetworKit::none;
}

NetworKit::index RooijBodlaenderSetCover::find_cardinality_frequency_set() {
    for (NetworKit::index i = 0; i < family.size(); i++) {
        auto &candidate = family.at(i);
        if (candidate.size() != 2) {
            continue;
        }
        bool satifies = std::all_of(
            candidate.begin(), candidate.end(),
            [this] (const auto& e) { return occurences.at(e).size() == 2; });
        if (satifies) {
            return i;
        }
    }
    return NetworKit::none;
}

}  /* namespace Koala */
