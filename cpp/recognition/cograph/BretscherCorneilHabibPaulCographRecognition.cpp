/*
 * BretscherCorneilHabibPaulCographRecognition.cpp
 *
 *  Created on: 21.04.2024
 *      Author: fixikmila
 */

#include <list>
#include <utility>
#include <vector>

#include "recognition/CographRecognition.hpp"

namespace Koala {

void BretscherCorneilHabibPaulCographRecognition::run() {
    hasRun = true;
    if (graph.numberOfNodes() <= 1) {
        is_cograph = State::COGRAPH;
        return;
    }
    std::vector<NetworKit::node> start(graph.nodeRange().begin(), graph.nodeRange().end());
    auto first = lex_bfs_minus(false, start);
    auto second = lex_bfs_minus(true, first.order);
    auto third = lex_bfs_minus(false, second.order);
    is_cograph = (neighbourhood_subset_property(true, second.order, second.slice_starts)
        && neighbourhood_subset_property(false, third.order, third.slice_starts))
        ? State::COGRAPH : State::NOT_COGRAPH;
}

BretscherCorneilHabibPaulCographRecognition::Info
BretscherCorneilHabibPaulCographRecognition::lex_bfs_minus(
        bool is_complement, const std::vector<NetworKit::node> &order) {
    std::vector<NetworKit::node> ans;
    std::vector<bool> used(graph.upperNodeIdBound());
    std::vector<std::vector<std::pair<NetworKit::node, NetworKit::count>>> what_ends_here(
        graph.upperNodeIdBound());
    std::vector<std::vector<NetworKit::count>> slice_starts(graph.upperNodeIdBound());
    std::list<std::list<NetworKit::node>> L;
    std::list<NetworKit::node> first;
    for (auto i : order) {
        first.push_back(i);
    }
    L.push_back(first);
    std::vector<std::list<std::list<NetworKit::node>>::iterator>
    in_which_list(graph.upperNodeIdBound(), L.begin()),
    previous_list(graph.upperNodeIdBound(), L.begin());
    std::vector<std::list<NetworKit::node>::iterator> in_which_position(
        graph.upperNodeIdBound());
    std::vector<bool> used_at_this_step(graph.upperNodeIdBound());
    auto it = L.front().begin();
    for (auto i : order) {
        in_which_position[i] = it++;
    }
    NetworKit::count i = 0;
    while (!L.empty()) {
        for (auto &[x, sz] : what_ends_here[i]) {
            NetworKit::count current_sum = 0;
            for (const auto &l : L) {
                if (l.size() + current_sum <= sz) {
                    slice_starts[x].emplace_back(i + current_sum);
                    current_sum += l.size();
                } else {
                    break;
                }
            }
        }

        auto x = L.front().front();
        auto slice_size = L.front().size();
        L.front().pop_front();
        if (L.front().empty()) {
            L.erase(L.begin());
        }
        ans.push_back(x);
        i++;
        used[x] = true;
        unsigned int SA = 0;
        for (auto v : graph.neighborRange(x)) {
            used_at_this_step[v] = true;
            previous_list[v] = in_which_list[v];
            if (used[v]) {
                continue;
            }
            auto number = in_which_list[v];
            if (slice_size != 1 && number == in_which_list[x]) {
                SA++;
            }
            auto previous = number;
            auto copy_number = number;
            auto end_it = L.end();
            end_it--;
            if (is_complement) {
                copy_number++;
            }
            if ((number == L.begin() && !is_complement) || (number == end_it && is_complement)) {
                L.insert(copy_number, std::list<NetworKit::node>());
            }
            if (!is_complement) {
                previous--;
            } else {
                previous++;
            }
            if (!previous->empty() && (!used_at_this_step[previous->front()]
                    || previous_list[previous->front()] != previous_list[v])) {
                std::list<NetworKit::node> insert;
                L.insert(copy_number, insert);
                previous = number;
                if (is_complement) {
                    previous++;
                } else {
                    previous--;
                }
            }
            previous->push_back(v);
            in_which_list[v] = previous;
            number->erase(in_which_position[v]);
            if (number->empty()) {
                L.erase(number);
            }
            auto end = previous->end();
            end--;
            in_which_position[v] = end;
        }
        if (is_complement) {
            SA = slice_size - 1 - SA;
        }
        for (auto v : graph.neighborRange(x)) {
            used_at_this_step[v] = false;
        }
        if (slice_size - SA - 1 != 0) {
            what_ends_here[i + SA].emplace_back(x, slice_size - SA - 1);
        }
    }
    return {slice_starts, ans};
}

bool BretscherCorneilHabibPaulCographRecognition::neighbourhood_subset_property(
        bool is_complement, const std::vector<NetworKit::node> &order,
        const std::vector<std::vector<NetworKit::count>> &slice_starts) {
    std::vector<bool> used(graph.upperNodeIdBound());
    std::vector<NetworKit::count> positions(graph.upperNodeIdBound(), NetworKit::none);
    for (NetworKit::count i = 0; i < order.size(); i++) {
        positions[order[i]] = i;
    }
    for (auto x : order) {
        const auto &starts = slice_starts[x];
        for (NetworKit::count j = 0; j + 1 < starts.size(); j++) {
            auto l = starts[j];
            auto r = starts[j + 1] - 1;
            auto y = order[l];
            auto p = order[starts[j + 1]];
            if (is_complement) {
                NetworKit::count count = 0;
                for (auto z : graph.neighborRange(p)) {
                    if (positions[z] < positions[p]) {
                        used[z] = true;
                    }
                    if (positions[z] >= l && positions[z] <= r) {
                        count++;
                    }
                }
                if (count != r - l + 1) {
                    return false;
                }
                for (auto z : graph.neighborRange(y)) {
                    if (positions[z] < positions[y] && !used[z]) {
                        return false;
                    }
                }
                for (auto z : graph.neighborRange(p)) {
                    if (positions[z] < positions[p]) {
                        used[z] = false;
                    }
                }
            } else {
                for (auto z : graph.neighborRange(y)) {
                    if (positions[z] < positions[y]) {
                        used[z] = true;
                    }
                }
                for (auto z : graph.neighborRange(p)) {
                    if (positions[z] < positions[p] && !used[z]) {
                        return false;
                    }
                }
                for (auto z : graph.neighborRange(y)) {
                    if (positions[z] < positions[y]) {
                        used[z] = false;
                    }
                }
            }
        }
    }
    return true;
}

}  // namespace Koala
