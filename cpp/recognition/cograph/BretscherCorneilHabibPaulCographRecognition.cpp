/*
 * BretscherCorneilHabibPaulCographRecognition.cpp
 *
 *  Created on: 2024
 *      Author: fixikmila
 */
#include <graph/GraphTools.hpp>

#include "recognition/CographRecognition.hpp"
#include "recognition/CoTree.hpp"

#include <list>

namespace Koala {

    bool BretscherCorneilHabibPaulCographRecognition::isCograph() const {
        assureFinished();
        return is_cograph;
    }

    void BretscherCorneilHabibPaulCographRecognition::run() {
        hasRun = true;
        std::vector<NetworKit::node> start;
        for (auto u : graph.nodeRange()) {
            start.push_back(u);
        }
        auto [bb, a] = LexBfsMinus(false, start);
        auto [borders1, b] = LexBfsMinus(true, a);  // graph complement
        auto [borders2, c] = LexBfsMinus(false, b);
        is_cograph = NeighbourhoodSubsetProperty(true, b, borders1) &&
            NeighbourhoodSubsetProperty(false, c, borders2);
    }

    std::pair<std::vector<std::vector<std::pair<int, int>>>, std::vector<NetworKit::node>>
        BretscherCorneilHabibPaulCographRecognition::
        LexBfsMinus(bool is_complement, std::vector<NetworKit::node> &a) {
        std::vector<NetworKit::node> ans;
        std::vector<bool> used(a.size());
        std::vector<std::vector<std::pair<NetworKit::node, unsigned int>>> what_ends_here(a.size());
        std::vector<std::vector<std::pair<int, int>>> borders(a.size());
        std::list<std::list<NetworKit::node>> L;
        std::list<NetworKit::node> first;
        for (auto i : a) {
            first.push_back(i);
        }
        L.push_back(first);
        std::vector<std::_List_iterator<std::list<NetworKit::node>>>
        in_which_list(a.size(), L.begin()), previous_list(a.size(), L.begin());
        std::vector<std::list<NetworKit::node>::iterator> in_which_position(a.size());
        std::vector<bool> used_at_this_step(a.size());
        auto it = L.front().begin();
        for (unsigned int i : a) {
            in_which_position[i] = it++;
        }
        int i = 0;
        while (!L.empty()) {
            for (auto &[x, sz] : what_ends_here[i]) {
                unsigned int current_sum = 0;
                for (const auto &l : L) {
                    if (l.size() + current_sum <= sz) {
                        borders[x].emplace_back(i + current_sum, i + current_sum + l.size() - 1);
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
                if ((number == L.begin() && !is_complement) ||
                (number == end_it && is_complement)) {
                    std::list<NetworKit::node> insert;
                    L.insert(copy_number, insert);
                }
                if (!is_complement) {
                    previous--;
                } else {
                    previous++;
                }
                if (!previous->empty() && (!graph.hasEdge(previous->front(), x)
                    || !used_at_this_step[previous->front()] ||
                    previous_list[previous->front()] != previous_list[v])) {
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
        return {borders, ans};
    }

    bool BretscherCorneilHabibPaulCographRecognition::NeighbourhoodSubsetProperty(
        bool is_complement, std::vector<NetworKit::node> a,
        std::vector<std::vector<std::pair<int, int>>> borders) {
        std::vector<bool> used(a.size());
        std::vector<int> positions(a.size());
        for (unsigned int i = 0; i < a.size(); i++) {
            positions[a[i]] = i;
        }
        for (auto i : a) {
            if (borders[i].empty()) {
                continue;
            }
            for (unsigned int j = 0; j < borders[i].size() - 1; j++) {
                auto [l, r] = borders[i][j];
                auto y = a[l];
                auto p = a[borders[i][j + 1].first];
                if (is_complement) {
                    int count = 0;
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
