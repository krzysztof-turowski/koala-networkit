/*
 * PlanarTree.hpp
 *
 *  Created on: 18.01.2021
 *      Author: Mikołaj Twaróg
 *      Ported by: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

template <typename Problem>
class PlanarTree {
 public:
    PlanarTree<Problem> *enclosing_tree;
    std::vector<Problem> t;
    int root, outer_face, enclosing_face, level;

 public:
    explicit PlanarTree(int size = 0, int l = 0) : level(l), enclosing_tree(nullptr) {
        for (int i = 0; i < size; i++) {
            t.emplace_back(l, this);
        }
    }

    PlanarTree(const PlanarTree<Problem> &other) {
        enclosing_tree = other.enclosing_tree, t = other.t;
        root = other.root, enclosing_face = other.enclosing_face, level = other.level;
        for (auto &p : t) {
            p.tree = this;
            if (!p.component_tree.empty()) {
                p.component_tree.enclosing_tree = this;
            }
        }
    }

    PlanarTree<Problem>& operator=(const PlanarTree<Problem> &other) {
        enclosing_tree = other.enclosing_tree, t = other.t;
        root = other.root, enclosing_face = other.enclosing_face, level = other.level;
        for (auto &p : t) {
            p.tree = this;
            if (!p.component_tree.empty()) {
                p.component_tree.enclosing_tree = this;
            }
        }
        return *this;
    }

    int merge(PlanarTree<Problem> &other, std::map<int, int> &place_in_comp, int node) {
        int s = size();
        for (int i = 0; i < other.size(); i++) {
            t.push_back(other[i]);
            t.back().tree = this;
            for (auto &c : t.back().children) {
                c += s;
            }
            t.back().parent += s;
        }

        int target = 0, other_v_root = other[other.root].label.first;
        std::queue<int> Q({node});
        while (!Q.empty()) {
            int i = Q.front();
            Q.pop();
            if (std::find(t[i].face.begin(), t[i].face.end(), other_v_root) != t[i].face.end()) {
                target = i;
                break;
            }
            for (int c : t[i].children) {
                Q.push(c);
            }
        }

        auto other_element = std::find_if(
            other[other.root].children.begin(), other[other.root].children.end(),
            [&] (auto child) { return other[child].label.second != other_v_root; });
        int other_place = place_in_comp[other[*other_element].label.second];

        int child = 0;
        for (int i = 0; i < t[target].children.size(); i++) {
            if (t[t[target].children[i]].label.second == other_v_root) {
                int t_place = INT16_MAX;
                if (i < t[target].children.size() - 1 && t[t[target].children[i + 1]].label.first == t[t[target].children[i + 1]].label.second) {
                    for (int child : t[t[target].children[i + 1]].children) {
                        if (t[child].label.second != other_v_root) {
                            t_place = place_in_comp[t[child].label.second];
                            break;
                        }
                    }
                }
                if (other_place < t_place) {
                    child = i + 1;
                    break;
                }
            }
            else if (t[t[target].children[i]].label.first == other_v_root) {
                child = i;
                break;
            }
        }
        auto t_element = std::find_if(
            t[target].children.begin(), t[target].children.end(),
            [&] (auto c) { return t[c].label.second != t[target].label.first; });
        int place = place_in_comp[t[*t_element].label.second];
        int shift = s + other.root;
        if (t[target].label.first == t[target].label.second && child == 0 && place < other_place) {
            t[target].children.push_back(shift);
        } else {
            t[target].children.insert(t[target].children.begin() + child, shift);
        }
        t[shift].parent = target;
        return shift;
    }

    Problem& operator[] (int x) {
        return t[x];
    }

    void remove_outer_face() {
        for (auto &p : t) {
            for (int &child : p.children) {
                if (child > outer_face) {
                    child--;
                }
            }
        }
        t.erase(t.begin() + outer_face);
        if (root > outer_face) {
            root--;
        }
    }

    Problem& get_enclosing_face() {
        return enclosing_tree->t[enclosing_face];
    }

    int size() {
        return t.size();
    }

    bool empty() {
        return t.empty();
    }

    void emplace_back() {
        t.emplace_back(level, this);
    }
};
