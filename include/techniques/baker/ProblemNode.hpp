/*
 * ProblemNode.hpp
 *
 *  Created on: 10.09.2023
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

#include <techniques/baker/PlanarTree.hpp>

struct Node {
    int parent, left, right;
    std::pair<int, int> label;
    std::vector<int> children, face;

    Node(): parent(-1) { }
};

class MaximizationNode : public Node {
 public:
    static int getInitial() {
        return 0;
    }

    static int getLimit(int k) {
        return k;
    }

    static int getNumberOfGraphs(int k, int max_level) {
        return ((max_level - 1) / (k + 1)) + 1;
    }

    static int getBest(int a, int b) {
        return std::max(a, b);
    }
    
    static void getLevels(
            const std::vector<int> &vertex_level, int r, int k, std::vector<std::vector<int>> &levels) {
        for (int i = 0; i < vertex_level.size(); i++) {
            int level = vertex_level[i] - 1;
            if (level % (k + 1) == r) {
                continue;
            }
            if (level < r) {
                levels[0].push_back(i);
            } else {
                levels[(level - r) / (k + 1)].push_back(i);
            }
        }
    }

};

class IndependentSetNode : public MaximizationNode {
 public:
    int level;
    std::vector<std::vector<int>> val;
    PlanarTree<IndependentSetNode> *tree, component_tree;

    IndependentSetNode(int l, PlanarTree<IndependentSetNode> *t = nullptr)
            : level(l), val(1 << l, std::vector<int>(1 << l)), tree(t) {
        if (l == 1) {
            val[0][0] = 0, val[0][1] = 1, val[1][0] = 1, val[1][1] = -INT16_MAX;
        }
    }

    IndependentSetNode(const IndependentSetNode &other) {
        parent = other.parent, left = other.left, right = other.right;
        label = other.label, children = other.children, face = other.face;
        level = other.level, val = other.val, tree = other.tree, component_tree = other.component_tree;
    }

    IndependentSetNode& operator=(const IndependentSetNode &other) {
        parent = other.parent, left = other.left, right = other.right;
        label = other.label, children = other.children, face = other.face;
        level = other.level, val = other.val, tree = other.tree, component_tree = other.component_tree;
        return *this;
    }

    IndependentSetNode& get_child(int x) {
        return (*tree)[children[x]];
    }

    void get_left_boundary(std::vector<int> &L) {
        L.push_back(label.first);
        if (tree->enclosing_tree == nullptr) {
            return;
        }
        auto &children = tree->get_enclosing_face().children;
        int child = left;
        if (left >= children.size()) {
            child--;
        }
        tree->enclosing_tree->t[children[child]].get_left_boundary(L);
    }

    void get_right_boundary(std::vector<int> &R) {
        R.push_back(label.second);
        if (tree->enclosing_tree == nullptr) {
            return;
        }
        tree->enclosing_tree->t[tree->get_enclosing_face().children[right - 1]].get_right_boundary(R);
    }

    void merge(std::vector<std::vector<int>> one, std::vector<std::vector<int>> other) {
        int count = 1 << level;
        for (int u = 0; u < count; u++){
            for (int v = 0; v < count; v++) {
                val[u][v] = -INT16_MAX;
                for (int z = 0; z < count; z++) {
                    int ones = 0;
                    for (int i = 0; i < level; i++) {
                        ones += (z & (1 << i)) > 0;
                    }
                    val[u][v] = std::max(val[u][v], one[u][z] + other[z][v] - ones);
                }
            }
        }
    }

    template<typename Graph>
    void adjust(Graph& g, std::set<std::pair<int, int>>& ae) {
        int count = 1 << (level - 1);
        if(label.first == label.second) {
            for (int u = 0; u < count; u++) {
                for (int v = 0; v < count; v++) {
                    val[(u << 1) + 1][(v << 1) + 1]--;
                    val[u << 1][(v << 1) + 1] = -INT16_MAX;
                    val[(u << 1) + 1][v << 1] = -INT16_MAX;
                }
            }
        }
        if (check_for_edge(label.first, label.second, g, ae)) {
            for (int u = 0; u < count; u++) {
                for (int v = 0; v < count; v++) {
                    val[(u << 1) + 1][(v << 1) + 1] = -INT16_MAX;
                }
            }
        }
    }

    void contract(IndependentSetNode& other, Graph& g, std::set<std::pair<int, int>>& ae) {
        int count = 1 << level;
        for (int u = 0; u < count; u++) {
            for (int v = 0; v < count; v++) {
                val[u][v] = std::max(other.val[(u << 1) + 1][(v << 1) + 1], other.val[u << 1][v << 1]);
            }
        }
    }

    template<typename Graph>
    IndependentSetNode extend(int z, Graph& g, std::set<std::pair<int, int>>& ae) {
        int count = 1 << level;
        IndependentSetNode res(level + 1);
        std::vector<int> L, R;
        get_left_boundary(L);
        get_right_boundary(R);
        for (int u = 0; u < count; u++) {
            for (int v = 0; v < count; v++) {
                res.val[u << 1][v << 1] = val[u][v];
                res.val[u << 1][(v << 1) + 1]= -INT16_MAX;
                res.val[(u << 1) + 1][v << 1] = -INT16_MAX;
                if (((u & 1) == 1 && check_for_edge(L[0], z, g, ae))
                    || ((v & 1) == 1 && check_for_edge(R[0], z, g, ae))){
                    res.val[(u << 1) + 1][(v << 1) + 1] = -INT16_MAX;
                } else {
                    res.val[(u << 1) + 1][(v << 1) + 1] = val[u][v] + 1;
                }
            }
        }

        return res;
    }

    template<typename Graph>
    void create(int child_num, Graph& g, std::set<std::pair<int, int>>& ae) {
        const std::vector<int>& children = tree->get_enclosing_face().children;
        std::vector<int> vertices;
        if (child_num < children.size()) {
            IndependentSetNode& child = tree->enclosing_tree->t[children[child_num]];
            child.get_left_boundary(vertices);
        } else {
            IndependentSetNode& child = tree->enclosing_tree->t[children[child_num - 1]];
            child.get_right_boundary(vertices);
        }

        int count = 1 << vertices.size();
        for (int u = 0; u < (count << 1); u++) {
            for (int v = 0; v < (count << 1); v++) {
                val[u][v] = -INT16_MAX;
            }
        }
        for (int i = 0; i < count; i++) {
            bool bad = false;
            for (int v = 0; v < vertices.size() - 1; v++) {
                if (((i >> v) & 1) && ((i >> (v + 1)) & 1) && check_for_edge(vertices[v], vertices[v+1], g, ae)) {
                    bad = true;
                    break;
                }
            }
            if (bad) {
                continue;
            }
            int ones = 0;
            for (int j = 0; j < vertices.size(); j++) {
                ones += (i & (1 << j)) > 0;
            }
            val[i << 1][i << 1] = ones;
            if ((i & 1) == 0 || !check_for_edge(vertices[0], label.first, g,  ae)) {
                val[(i << 1) + 1][i << 1] = ones + 1;
            }
            if ((i & 1) == 0 || !check_for_edge(vertices[0], label.second, g, ae)) {
                val[i << 1][(i << 1) + 1] = ones + 1;
            }
            if ((i & 1) == 0 || (!check_for_edge(vertices[0], label.second, g, ae)
                    && !check_for_edge(vertices[0], label.first, g,  ae))) {
                val[(i << 1) + 1][(i << 1) + 1] = ones + 2;
            }
        }
    }

    int result() {
        return std::max(val[0][0], val[1][1]);
    }
};

class MinimizationNode : public Node {
 public:
    static int getInitial() {
        return INT16_MAX;
    }

    static int getLimit(int k) {
        return k + 1;
    }

    static int getNumberOfGraphs(int k, int max_level) {
        return ((max_level - 1) / k) + 2;
    }

    static int getBest(int a, int b) {
        return std::min(a, b);
    }
    
    static void getLevels(
            const std::vector<int> &vertex_level, int r, int k, std::vector<std::vector<int>> &levels) {
        for (int i = 0; i < vertex_level.size(); i++) {
            int level = vertex_level[i] - 1;
            if (level < r) {
                levels[0].push_back(i);
                continue;
            }
            if (level % k == r) {
                levels[(level - r) / k].push_back(i);
            }
            levels[((level - r) / k) + 1].push_back(i);
        }
    }
};

class VertexCoverNode : public MinimizationNode {
 public:
    int level;
    std::vector<std::vector<int>> val;
    PlanarTree<VertexCoverNode> *tree, component_tree;

    VertexCoverNode(int l): VertexCoverNode(l, nullptr) {}

    VertexCoverNode(int l, PlanarTree<VertexCoverNode> *t = nullptr):
            level(l), val(1 << l,std::vector<int>(1 << l)), tree(t) {
        if(l == 1) {
            val[0][0] = INT16_MAX - 1;
            val[0][1] = val[1][0] = 1;
            val[1][1] = 2;
        }
    }

    VertexCoverNode(const VertexCoverNode& other) {
        parent = other.parent, left = other.left, right = other.right;
        label = other.label, children = other.children, face = other.face;
        level = other.level, val = other.val, tree = other.tree, component_tree = other.component_tree;
    }

    VertexCoverNode& operator=(const VertexCoverNode& other) {
        parent = other.parent, left = other.left, right = other.right;
        label = other.label, children = other.children, face = other.face;
        level = other.level, val = other.val, tree = other.tree, component_tree = other.component_tree;
        return *this;
    }

    VertexCoverNode& get_child(int x) {
        return (*tree)[children[x]];
    }

    void get_left_boundary(std::vector<int> &L) {
        L.push_back(label.first);
        if (tree->enclosing_tree == nullptr) {
            return;
        }
        auto &children = tree->get_enclosing_face().children;
        int child = left;
        if (left >= children.size()) {
            child--;
        }
        tree->enclosing_tree->t[children[child]].get_left_boundary(L);
    }

    void get_right_boundary(std::vector<int> &R) {
        R.push_back(label.second);
        if (tree->enclosing_tree == nullptr) {
            return;
        }
        tree->enclosing_tree->t[tree->get_enclosing_face().children[right - 1]].get_right_boundary(R);
    }

    void merge(std::vector<std::vector<int>> one, std::vector<std::vector<int>> other) {
        int count = 1 << level;
        for (int u = 0; u < count; u++){
            for (int v = 0; v < count; v++) {
                val[u][v] = INT16_MAX - 1;
                for (int z = 0; z < count; z++) {
                    int ones = 0;
                    for (int i = 0; i < level; i++) {
                        ones += (z & (1 << i)) > 0;
                    }

                    val[u][v] = std::min(val[u][v], one[u][z] + other[z][v] - ones);
                }
            }
        }
    }

    template<typename Graph>
    void adjust(Graph& g, std::set<std::pair<int, int>>& ae) {
        int count = 1 << (level - 1);

        if(label.first == label.second) {
            for (int u = 0; u < count; u++) {
                for (int v = 0; v < count; v++) {
                    val[(u << 1) + 1][(v << 1) + 1]--;
                    val[u << 1][(v << 1) + 1] = INT16_MAX - 1;
                    val[(u << 1) + 1][v << 1] = INT16_MAX - 1;
                }
            }
        }
        if (check_for_edge(label.first, label.second, g, ae)) {
            for (int u = 0; u < count; u++) {
                for (int v = 0; v < count; v++) {
                    val[u << 1][v << 1] = INT16_MAX - 1;
                }
            }
        }
    }

    void contract(VertexCoverNode& other, Graph& g, std::set<std::pair<int, int>>& ae) {
        int count = 1 << level;
        for (int u = 0; u < count; u++) {
            for (int v = 0; v < count; v++) {
                val[u][v] = std::min(other.val[(u << 1) + 1][(v << 1) + 1], other.val[u << 1][v << 1]);
            }
        }
    }

    template<typename Graph>
    VertexCoverNode extend(int z, Graph& g, std::set<std::pair<int, int>>& ae) {
        int count = 1 << level;
        VertexCoverNode res(level + 1);
        std::vector<int> L, R;
        get_left_boundary(L), get_right_boundary(R);
        for (int u = 0; u < count; u++) {
            for (int v = 0; v < count; v++) {
                res.val[(u << 1) + 1][(v << 1) + 1] = val[u][v] + 1;
                res.val[u << 1][(v << 1) + 1]= INT16_MAX - 1;
                res.val[(u << 1) + 1][v << 1] = INT16_MAX - 1;

                if (((u & 1) == 1 || !check_for_edge(L[0], z, g, ae))
                    && ((v & 1) == 1 || !check_for_edge(R[0], z, g, ae))){
                    res.val[u << 1][v << 1] = val[u][v];
                } else {
                    res.val[u << 1][v << 1] = INT16_MAX - 1;
                }
            }
        }

        return res;
    }

    template<typename Graph>
    void create(int child_num, Graph& g, std::set<std::pair<int, int>>& ae) {
        const std::vector<int>& children = tree->get_enclosing_face().children;
        std::vector<int> vertices;

        if (child_num < children.size()) {
            VertexCoverNode& child = tree->enclosing_tree->t[children[child_num]];
            child.get_left_boundary(vertices);
        } else {
            VertexCoverNode& child = tree->enclosing_tree->t[children[child_num - 1]];
            child.get_right_boundary(vertices);
        }

        int count = 1 << vertices.size();

        for (int u = 0; u < (count << 1); u++) {
            for (int v = 0; v < (count << 1); v++) {
                val[u][v] = INT16_MAX - 1;
            }
        }

        for (int i = 0; i < count; i++) {

            bool bad = false;
            for (int v = 0; v < vertices.size() - 1; v++) {
                if (!((i >> v) & 1) && !((i >> (v + 1)) & 1) && check_for_edge(vertices[v], vertices[v+1], g, ae)) {
                    bad = true;
                    break;
                }
            }

            if (bad) {
                continue;
            }

            int ones = 0;
            for (int j = 0; j < vertices.size(); j++) {
                ones += (i & (1 << j)) > 0;
            }

            val[(i << 1) + 1][(i << 1) + 1] = ones + 2;

            if ((i & 1) == 1 || !check_for_edge(vertices[0], label.first, g,  ae)) {
                val[i << 1][(i << 1) + 1] = ones + 1;
            }

            if ((i & 1) == 1 || !check_for_edge(vertices[0], label.second, g, ae)) {
                val[(i << 1) + 1][i << 1] = ones + 1;
            }

            if ((i & 1) == 1 || (!check_for_edge(vertices[0], label.second, g, ae)
                                 && !check_for_edge(vertices[0], label.first, g,  ae))) {
                val[i << 1][i << 1] = ones;
            }
        }
    }

    int result() {
        return std::min(val[0][0], val[1][1]);
    }
};

struct DominatingSetNode : public MinimizationNode {
 public:
    int level;
    std::vector<std::vector<int>> val;
    PlanarTree<DominatingSetNode> *tree, component_tree;

    DominatingSetNode(int l, PlanarTree<DominatingSetNode> *t = nullptr)
            : level(l), val(pow(3, l), std::vector<int>(pow(3, l))), tree(t) {
        if(l == 1) {
            val[0][0] = INT16_MAX - 1;
            val[0][1] = 1;
            val[0][2] = INT16_MAX - 1;
            val[1][0] = 1;
            val[1][1] = 2;
            val[1][2] = 1;
            val[2][0] = INT16_MAX - 1;
            val[2][1] = 1;
            val[2][2] = 0;
        }
    }

    DominatingSetNode(const DominatingSetNode& other) {
        parent = other.parent, left = other.left, right = other.right;
        label = other.label, children = other.children, face = other.face;
        level = other.level, val = other.val, tree = other.tree, component_tree = other.component_tree;
    }

    DominatingSetNode& operator=(const DominatingSetNode& other) {
        parent = other.parent, left = other.left, right = other.right;
        label = other.label, children = other.children, face = other.face;
        level = other.level, val = other.val, tree = other.tree, component_tree = other.component_tree;
        return *this;
    }

    DominatingSetNode& get_child(int x) {
        return (*tree)[children[x]];
    }

    void get_left_boundary(std::vector<int> &L) {
        L.push_back(label.first);
        if (tree->enclosing_tree == nullptr) {
            return;
        }
        auto &children = tree->get_enclosing_face().children;
        int child = left;
        if (left >= children.size()) {
            child--;
        }
        tree->enclosing_tree->t[children[child]].get_left_boundary(L);
    }

    void get_right_boundary(std::vector<int> &R) {
        R.push_back(label.second);
        if (tree->enclosing_tree == nullptr) {
            return;
        }
        tree->enclosing_tree->t[tree->get_enclosing_face().children[right - 1]].get_right_boundary(R);
    }

    void merge(std::vector<std::vector<int>> one, std::vector<std::vector<int>> other) {
        int count = val.size();
        for (int u = 0; u < count; u++){
            for (int v = 0; v < count; v++) {
                val[u][v] = INT16_MAX - 1;
                for (int one_z = 0; one_z < count; one_z++) {
                    int other_z = one_z, three_z = one_z, ones = 0;
                    for (int i = one_z, j = 2, it = 0; it < level; i /= 3, j *= 3, it++) {
                        int mod = i % 3;
                        if (mod == 0) {
                            other_z += j;
                        } else if (mod == 1) {
                            ones++;
                        } else if (mod == 2) {
                            other_z -= j, three_z -= j;
                        }
                    }
                    val[u][v] = std::min(val[u][v], one[u][one_z] + other[three_z][v] - ones);
                    val[u][v] = std::min(val[u][v], one[u][three_z] + other[one_z][v] - ones);
                    val[u][v] = std::min(val[u][v], one[u][one_z] + other[other_z][v] - ones);
                }
            }
        }
    }

    template<typename Graph>
    void adjust(Graph& g, std::set<std::pair<int, int>>& ae) {
        int count = pow(3, level - 1);
        if(label.first == label.second) {
            for (int u = 0; u < count; u++) {
                for (int v = 0; v < count; v++) {
                    val[u * 3][v * 3] = std::min(val[u * 3][v * 3], val[u * 3][(v * 3) + 2]);
                    val[u * 3][v * 3] = std::min(val[u * 3][v * 3], val[(u * 3) + 2][v * 3]);
                    val[(u * 3) + 2][(v * 3) + 2] = std::min(val[(u * 3) + 2][(v * 3) + 2], val[u * 3][v * 3]);
                    val[(u * 3) + 1][(v * 3) + 1]--;
                    val[u * 3][(v * 3) + 1] = INT16_MAX - 1;
                    val[u * 3][(v * 3) + 2] = val[(u * 3) + 2][(v * 3) + 2];
                    val[(u * 3) + 1][v * 3] = INT16_MAX - 1;
                    val[(u * 3) + 1][(v * 3) + 2] = INT16_MAX - 1;
                    val[(u * 3) + 2][v * 3] = val[(u * 3) + 2][(v * 3) + 2];
                    val[(u * 3) + 2][(v * 3) + 1] = INT16_MAX - 1;
                }
            }
        }
        if (check_for_edge(label.first, label.second, g, ae)) {
            for (int u = 0; u < count; u++) {
                for (int v = 0; v < count; v++) {
                    val[(u * 3) + 1][v * 3] = std::min(val[(u * 3) + 1][v * 3], val[(u * 3) + 1][(v * 3) + 2]);
                    val[u * 3][(v * 3) + 1] = std::min(val[u * 3][(v * 3) + 1], val[(u * 3) + 2][(v * 3) + 1]);
                }
            }
        }
    }

    void contract(DominatingSetNode& other, Graph& g, std::set<std::pair<int, int>>& ae) {
        int count = val.size();
        for (int u = 0; u < count; u++) {
            for (int v = 0; v < count; v++) {
                val[u][v] = std::min(other.val[(u * 3) + 1][(v * 3) + 1], other.val[u * 3][v * 3]);
                if ((u % 3 == 1 && check_for_edge(label.first, other.label.first, g, ae)) ||
                (v % 3 == 1 && check_for_edge(label.second, other.label.first, g, ae))) {
                    val[u][v] = std::min(val[u][v], other.val[(u * 3) + 2][(v * 3) + 2]);
                }
            }
        }
    }

    template<typename Graph>
    DominatingSetNode extend(int z, Graph& g, std::set<std::pair<int, int>>& ae) {
        int count = val.size();
        DominatingSetNode res(level + 1);
        std::vector<int> L, R;
        get_left_boundary(L);
        get_right_boundary(R);
        for (int u = 0; u < count; u++) {
            for (int v = 0; v < count; v++) {
                res.val[u * 3][(v * 3) + 1]= INT16_MAX - 1;
                res.val[(u * 3) + 1][v * 3] = INT16_MAX - 1;
                res.val[(u * 3) + 1][(v * 3) + 2] = INT16_MAX - 1;
                res.val[(u * 3) + 2][(v * 3) + 1] = INT16_MAX - 1;
                res.val[(u * 3) + 2][(v * 3) + 2] = val[u][v];
                if (((u % 3) == 1 && check_for_edge(L[0], z, g, ae))
                        || ((v % 3) == 1 && check_for_edge(R[0], z, g, ae))){
                    res.val[u * 3][v * 3] = val[u][v];
                } else {
                    res.val[u * 3][v * 3] = INT16_MAX - 1;
                }
                res.val[u * 3][(v * 3) + 2]= res.val[(u * 3) + 2][(v * 3) + 2];
                res.val[(u * 3) + 2][v * 3] = res.val[(u * 3) + 2][(v * 3) + 2];
                res.val[(u * 3) + 1][(v * 3) + 1] = val[u][v];
                int u_temp = u, v_temp = v;
                if ((u % 3) == 0 && check_for_edge(L[0], z, g, ae)) {
                    u_temp += 2;
                    res.val[(u * 3) + 1][(v * 3) + 1] = std::min(res.val[(u * 3) + 1][(v * 3) + 1], val[u_temp][v]);
                }
                if ((v % 3) == 0 && check_for_edge(R[0], z, g, ae)) {
                    v_temp += 2;
                    res.val[(u * 3) + 1][(v * 3) + 1] = std::min(res.val[(u * 3) + 1][(v * 3) + 1], val[u][v_temp]);
                }
                res.val[(u * 3) + 1][(v * 3) + 1] = std::min(res.val[(u * 3) + 1][(v * 3) + 1], val[u_temp][v_temp]);
                res.val[(u * 3) + 1][(v * 3) + 1]++;
            }
        }

        return res;
    }

    template<typename Graph>
    void create(int child_num, Graph& g, std::set<std::pair<int, int>>& ae) {
        const std::vector<int>& children = tree->get_enclosing_face().children;
        std::vector<int> vertices;

        if (child_num < children.size()) {
            DominatingSetNode& child = tree->enclosing_tree->t[children[child_num]];
            child.get_left_boundary(vertices);
        } else {
            DominatingSetNode& child = tree->enclosing_tree->t[children[child_num - 1]];
            child.get_right_boundary(vertices);
        }
        int count = pow(3, vertices.size());
        for (int u = 0; u < (count * 3); u++) {
            for (int v = 0; v < (count * 3); v++) {
                val[u][v] = INT16_MAX - 1;
            }
        }
        for (int i = 0; i < count; i++) {
            for (int j = 0; j < count; j++) {
                bool bad = false;
                int ones = 0;
                for (int v = 0, p = 1; v < vertices.size(); v++, p *= 3) {
                    if (((i / p) % 3 == 1) != ((j / p) % 3 == 1)) {
                        bad = true;
                        break;
                    }
                    ones += (i / p) % 3 == 1;
                }

                for (int v = 1, p = 3; v < vertices.size(); v++, p *= 3) {
                    bool up_dominated = false;
                    if (v < vertices.size() - 1) {
                        up_dominated = (i / (p * 3)) % 3 == 1 && check_for_edge(vertices[v], vertices[v + 1], g, ae);
                    }
                    bool down_dominated = (i / (p / 3)) % 3 == 1 && check_for_edge(vertices[v], vertices[v - 1], g, ae);
                    if ((i / p) % 3 == 0 && (j / p) % 3 == 0 &&
                        !((i / (p * 3)) % 3 == 1 && check_for_edge(vertices[v], vertices[v + 1], g, ae)
                          || (i / (p / 3)) % 3 == 1 && check_for_edge(vertices[v], vertices[v - 1], g, ae))) {
                        bad = true;
                        break;
                    }
                }
                if (bad) {
                    continue;
                }
                bool x_edge = check_for_edge(vertices[0], label.first, g, ae);
                bool y_edge = check_for_edge(vertices[0], label.second, g, ae);
                val[(i * 3) + 1][(j * 3) + 1] = ones + 2;
                val[(i * 3) + 2][(j * 3) + 2] = ones;
                val[(i * 3) + 1][(j * 3) + 2] = ones + 1;
                val[(i * 3) + 2][(j * 3) + 1] = ones + 1;
                if (x_edge && i % 3 == 1) {
                    val[i * 3][(j * 3) + 1] = ones + 1;
                    val[i * 3][(j * 3) + 2] = ones;
                }
                if (y_edge && i % 3 == 1) {
                    val[(i * 3) + 1][j * 3] = ones + 1;
                    val[(i * 3) + 2][j * 3] = ones;
                }
                if (x_edge && y_edge && i % 3 == 1) {
                    val[i * 3][j * 3] = ones;
                }
                bool z_edge = !(vertices.size() == 1 || !check_for_edge(vertices[0], vertices[1], g, ae));
                if (i % 3 == 0 && j % 3 == 0 && !((i/ 3) % 3 == 1 && z_edge)) {
                    val[(i * 3) + 2][(j * 3) + 2] = INT16_MAX;
                    if (!x_edge) {
                        val[(i * 3) + 1][(j * 3) + 2] = INT16_MAX;
                    }
                    if (!y_edge) {
                        val[(i * 3) + 2][(j * 3) + 1] = INT16_MAX;
                    }
                    if (!x_edge && !y_edge) {
                        val[(i * 3) + 1][(j * 3) + 1] = INT16_MAX;
                    }
                }
            }
        }
    }

    int result() {
        return std::min(val[0][0], val[1][1]);
    }
};
