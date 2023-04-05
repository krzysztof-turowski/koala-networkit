//
// Created by mikolaj on 1/18/21.
//

#ifndef TECHNIKA_BAKER_PROBLEMS_HPP
#define TECHNIKA_BAKER_PROBLEMS_HPP

//#include "../utils/level_face_traversal.hpp"

template <typename Problem>
struct tree{
    tree<Problem>* enclosing_tree;
    int enclosing_face;
    std::vector<Problem> t;
    int root;
    int outer_face;
    int level;

    tree():enclosing_tree(nullptr){}
    explicit tree(int size): t(size), enclosing_tree(nullptr){}
    explicit tree(int size, int l): level(l), enclosing_tree(nullptr){
        for (int i = 0; i < size; i++) {
            t.emplace_back(l, this);
//            t.back().my_tree = this;
        }
    }
    tree(const tree<Problem>& t2) {
        t = t2.t;
        enclosing_tree = t2.enclosing_tree;
        enclosing_face = t2.enclosing_face;
        root = t2.root;
        level = t2.level;
        for (auto& p : t) {
            p.my_tree = this;
            if (!p.component_tree.empty()) {
                p.component_tree.enclosing_tree = this;
            }
        }
    }

    Problem& operator[] (int x) {
        return t[x];
    }

    tree<Problem>& operator= (const tree<Problem>& t2) {
        t = t2.t;
        enclosing_tree = t2.enclosing_tree;
        enclosing_face = t2.enclosing_face;
        root = t2.root;
        level = t2.level;
        for (auto& p : t) {
            p.my_tree = this;
            if (!p.component_tree.empty()) {
                p.component_tree.enclosing_tree = this;
            }
        }
        return *this;
    }

    int merge(tree<Problem>& t2, std::map<int, int>& place_in_comp, int starting_node) {
        int s = size();
        for (int i = 0; i < t2.size(); i++) {
            t.push_back(t2[i]);
            t.back().my_tree = this;

            for (int& c : t.back().children) {
                c += s;
            }

            t.back().parent += s;
        }

        int target = 0;
        int t2_v_root = t2[t2.root].label.first;
        std::queue<int> q;
        q.push(starting_node);

        while (!q.empty()) {
            int i = q.front();
            q.pop();

            if (std::find(t[i].face.begin(), t[i].face.end(), t2_v_root) != t[i].face.end()) {
                target = i;
                break;
            }

            for (int c : t[i].children) {
                q.push(c);
            }
        }

        int t2_place;

        for (int child : t2[t2.root].children) {
            if (t2[child].label.second != t2_v_root) {
                t2_place = place_in_comp[t2[child].label.second];
                break;
            }
        }

        int child_num = 0;
        for (int i = 0; i < t[target].children.size(); i++) {
            if (t[t[target].children[i]].label.second == t2_v_root) {
                int t1_place = INT16_MAX;

                if (i < t[target].children.size() - 1 && t[t[target].children[i + 1]].label.first == t[t[target].children[i + 1]].label.second) {
                    for (int child : t[t[target].children[i + 1]].children) {
                        if (t[child].label.second != t2_v_root) {
                            t1_place = place_in_comp[t[child].label.second];
                            break;
                        }
                    }
                }

                if (t2_place < t1_place) {
                    child_num = i + 1;
                    break;
                }
            }
            else if (t[t[target].children[i]].label.first == t2_v_root) {
                child_num = i;
                break;
            }
        }

        int t1_place;

        for (int child : t[target].children) {
            if (t[child].label.second != t[target].label.first) {
                t1_place = place_in_comp[t[child].label.second];
                break;
            }
        }

        if (t[target].label.first == t[target].label.second && child_num == 0 && t1_place < t2_place) {
            t[target].children.push_back(s + t2.root);
        } else {
            t[target].children.insert(t[target].children.begin() + child_num, s + t2.root);
        }

        t[s + t2.root].parent = target;

        return s + t2.root;
    }

    void remove_outer_face() {
        for (auto& p : t) {
            for (int& child : p.children) {
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
//        t.back().my_tree = this;
    }
};

struct node
{
    int parent;
    std::pair<int, int> label;
    int LB;
    int RB;
    std::vector<int> children;
    std::vector<int> face;
    node(): parent(-1) {}
};

struct independent_set : node
{
    tree<independent_set>* my_tree;
    tree<independent_set> component_tree;
    std::vector< std::vector<int> > val;
    int level;

    independent_set(int l): independent_set(l, nullptr) {}

    independent_set(int l, tree<independent_set>* mt):
    level(l), val(1 << l,std::vector<int>(1 << l)), my_tree(mt) {
        if(l == 1) {
            val[0][0] = 0;
            val[0][1] = 1;
            val[1][0] = 1;
            val[1][1] = -INT16_MAX;
        }
    }

    independent_set(const independent_set& two): my_tree(two.my_tree) {
        parent = two.parent;
        label = two.label;
        LB = two.LB;
        RB = two.RB;
        children = two.children;
        face = two.face;
        val = two.val;
        component_tree = two.component_tree;
        level = two.level;
//        my_tree = two.my_tree;
    }

    independent_set& operator=(const independent_set& two) {
        parent = two.parent;
        label = two.label;
        LB = two.LB;
        RB = two.RB;
        children = two.children;
        face = two.face;
        val = two.val;
        component_tree = two.component_tree;
        my_tree = two.my_tree;
        level = two.level;
        return *this;
    }

    independent_set& get_child(int x) {
        return my_tree->t[children[x]];
    }

    // val = {v_0, ..., v_(1 << (level*2)}
    // v_0 <- x
    // v_(1 << level) <- y

    void get_left_boundary(std::vector<int>& lb) {
        lb.push_back(label.first);

        if (my_tree->enclosing_tree == nullptr)
            return;

        auto& children = my_tree->enclosing_tree->t[my_tree->enclosing_face].children;
        int child = LB;
        if (LB >= children.size()) {
            child--;
        }
        my_tree->enclosing_tree->t[children[child]].get_left_boundary(lb);
    }

    void get_right_boundary(std::vector<int>& rb) {
        rb.push_back(label.second);

        if (my_tree->enclosing_tree == nullptr)
            return;

        my_tree->enclosing_tree->t[my_tree->enclosing_tree->t[my_tree->enclosing_face].children[RB - 1]].get_right_boundary(rb);
    }

    void merge(std::vector< std::vector<int> > one, std::vector< std::vector<int> > two) {
//        std::vector<int> copy(val);

        int count = 1 << level;

        for (int u = 0; u < count; u++){
            for (int v = 0; v < count; v++) {
                val[u][v] = -INT16_MAX;
                for (int z = 0; z < count; z++) {
                    int ones = 0;
                    for (int i = 0; i < level; i++) {
                        ones += (z & (1 << i)) > 0;
                    }

                    val[u][v] = std::max(val[u][v], one[u][z] + two[z][v] - ones);
                }
            }
        }
    }

    template<typename Graph>
    void adjust(Graph& g, std::set<std::pair<int, int> >& ae) {
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

    void contract(independent_set& two, Graph& g, std::set< std::pair<int, int> >& ae) {
        int count = 1 << level;

        for (int u = 0; u < count; u++) {
            for (int v = 0; v < count; v++) {
                val[u][v] = std::max(two.val[(u << 1) + 1][(v << 1) + 1], two.val[u << 1][v << 1]);
            }
        }
    }

    template<typename Graph>
    independent_set extend(int z, Graph& g, std::set<std::pair<int, int> >& ae) {
        int count = 1 << level;

        independent_set res(level + 1);

        std::vector<int> lb;
        std::vector<int> rb;

        get_left_boundary(lb);
        get_right_boundary(rb);

        for (int u = 0; u < count; u++) {
            for (int v = 0; v < count; v++) {
                res.val[u << 1][v << 1] = val[u][v];
                res.val[u << 1][(v << 1) + 1]= -INT16_MAX;
                res.val[(u << 1) + 1][v << 1] = -INT16_MAX;

                if (((u & 1) == 1 && check_for_edge(lb[0], z, g, ae))
                    || ((v & 1) == 1 && check_for_edge(rb[0], z, g, ae))){
                    res.val[(u << 1) + 1][(v << 1) + 1] = -INT16_MAX;
                } else {
                    res.val[(u << 1) + 1][(v << 1) + 1] = val[u][v] + 1;
                }
            }
        }

        return res;
    }

    template<typename Graph>
    void create(int child_num, Graph& g, std::set<std::pair<int, int> >& ae) {
        const std::vector<int>& children = my_tree->enclosing_tree->t[my_tree->enclosing_face].children;
        std::vector<int> vertices;

        if (child_num < children.size()) {
            independent_set& child = my_tree->enclosing_tree->t[children[child_num]];
            child.get_left_boundary(vertices);
        } else {
            independent_set& child = my_tree->enclosing_tree->t[children[child_num - 1]];
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


struct vertex_cover : node
{
    tree<vertex_cover>* my_tree;
    tree<vertex_cover> component_tree;
    std::vector< std::vector<int> > val;
    int level;

    vertex_cover(int l): vertex_cover(l, nullptr) {}

    vertex_cover(int l, tree<vertex_cover>* mt):
            level(l), val(1 << l,std::vector<int>(1 << l)), my_tree(mt) {
        if(l == 1) {
            val[0][0] = INT16_MAX - 1;
            val[0][1] = 1;
            val[1][0] = 1;
            val[1][1] = 2;
        }
    }

    vertex_cover(const vertex_cover& two): my_tree(two.my_tree) {
        parent = two.parent;
        label = two.label;
        LB = two.LB;
        RB = two.RB;
        children = two.children;
        face = two.face;
        val = two.val;
        component_tree = two.component_tree;
        level = two.level;
//        my_tree = two.my_tree;
    }

    vertex_cover& operator=(const vertex_cover& two) {
        parent = two.parent;
        label = two.label;
        LB = two.LB;
        RB = two.RB;
        children = two.children;
        face = two.face;
        val = two.val;
        component_tree = two.component_tree;
        my_tree = two.my_tree;
        level = two.level;
        return *this;
    }

    vertex_cover& get_child(int x) {
        return my_tree->t[children[x]];
    }

    // val = {v_0, ..., v_(1 << (level*2)}
    // v_0 <- x
    // v_(1 << level) <- y

    void get_left_boundary(std::vector<int>& lb) {
        lb.push_back(label.first);

        if (my_tree->enclosing_tree == nullptr)
            return;

        auto& children = my_tree->enclosing_tree->t[my_tree->enclosing_face].children;
        int child = LB;
        if (LB >= children.size()) {
            child--;
        }
        my_tree->enclosing_tree->t[children[child]].get_left_boundary(lb);
    }

    void get_right_boundary(std::vector<int>& rb) {
        rb.push_back(label.second);

        if (my_tree->enclosing_tree == nullptr)
            return;

        my_tree->enclosing_tree->t[my_tree->enclosing_tree->t[my_tree->enclosing_face].children[RB - 1]].get_right_boundary(rb);
    }

    void merge(std::vector< std::vector<int> > one, std::vector< std::vector<int> > two) {
//        std::vector<int> copy(val);

        int count = 1 << level;

        for (int u = 0; u < count; u++){
            for (int v = 0; v < count; v++) {
                val[u][v] = INT16_MAX - 1;
                for (int z = 0; z < count; z++) {
                    int ones = 0;
                    for (int i = 0; i < level; i++) {
                        ones += (z & (1 << i)) > 0;
                    }

                    val[u][v] = std::min(val[u][v], one[u][z] + two[z][v] - ones);
                }
            }
        }
    }

    template<typename Graph>
    void adjust(Graph& g, std::set<std::pair<int, int> >& ae) {
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

    void contract(vertex_cover& two, Graph& g, std::set< std::pair<int, int> >& ae) {
        int count = 1 << level;

        for (int u = 0; u < count; u++) {
            for (int v = 0; v < count; v++) {
                val[u][v] = std::min(two.val[(u << 1) + 1][(v << 1) + 1], two.val[u << 1][v << 1]);
            }
        }
    }

    template<typename Graph>
    vertex_cover extend(int z, Graph& g, std::set< std::pair<int, int> >& ae) {
        int count = 1 << level;

        vertex_cover res(level + 1);

        std::vector<int> lb;
        std::vector<int> rb;

        get_left_boundary(lb);
        get_right_boundary(rb);

        for (int u = 0; u < count; u++) {
            for (int v = 0; v < count; v++) {
                res.val[(u << 1) + 1][(v << 1) + 1] = val[u][v] + 1;
                res.val[u << 1][(v << 1) + 1]= INT16_MAX - 1;
                res.val[(u << 1) + 1][v << 1] = INT16_MAX - 1;

                if (((u & 1) == 1 || !check_for_edge(lb[0], z, g, ae))
                    && ((v & 1) == 1 || !check_for_edge(rb[0], z, g, ae))){
                    res.val[u << 1][v << 1] = val[u][v];
                } else {
                    res.val[u << 1][v << 1] = INT16_MAX - 1;
                }
            }
        }

        return res;
    }

    template<typename Graph>
    void create(int child_num, Graph& g, std::set<std::pair<int, int> >& ae) {
        const std::vector<int>& children = my_tree->enclosing_tree->t[my_tree->enclosing_face].children;
        std::vector<int> vertices;

        if (child_num < children.size()) {
            vertex_cover& child = my_tree->enclosing_tree->t[children[child_num]];
            child.get_left_boundary(vertices);
        } else {
            vertex_cover& child = my_tree->enclosing_tree->t[children[child_num - 1]];
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


struct dominating_set : node
{
    tree<dominating_set>* my_tree;
    tree<dominating_set> component_tree;
    std::vector< std::vector<int> > val;
    int level;

    dominating_set(int l): dominating_set(l, nullptr) {}

    dominating_set(int l, tree<dominating_set>* mt):
            level(l), val(pow(3, l),std::vector<int>(pow(3, l))), my_tree(mt) {
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

    dominating_set(const dominating_set& two): my_tree(two.my_tree) {
        parent = two.parent;
        label = two.label;
        LB = two.LB;
        RB = two.RB;
        children = two.children;
        face = two.face;
        val = two.val;
        component_tree = two.component_tree;
        level = two.level;
//        my_tree = two.my_tree;
    }

    dominating_set& operator=(const dominating_set& two) {
        parent = two.parent;
        label = two.label;
        LB = two.LB;
        RB = two.RB;
        children = two.children;
        face = two.face;
        val = two.val;
        component_tree = two.component_tree;
        my_tree = two.my_tree;
        level = two.level;
        return *this;
    }

    dominating_set& get_child(int x) {
        return my_tree->t[children[x]];
    }

    // val = {v_0, ..., v_(1 << (level*2)}
    // v_0 <- x
    // v_(1 << level) <- y

    void get_left_boundary(std::vector<int>& lb) {
        lb.push_back(label.first);

        if (my_tree->enclosing_tree == nullptr)
            return;

        auto& children = my_tree->enclosing_tree->t[my_tree->enclosing_face].children;
        int child = LB;
        if (LB >= children.size()) {
            child--;
        }
        my_tree->enclosing_tree->t[children[child]].get_left_boundary(lb);
    }

    void get_right_boundary(std::vector<int>& rb) {
        rb.push_back(label.second);

        if (my_tree->enclosing_tree == nullptr)
            return;

        my_tree->enclosing_tree->t[my_tree->enclosing_tree->t[my_tree->enclosing_face].children[RB - 1]].get_right_boundary(rb);
    }

    void merge(std::vector< std::vector<int> > one, std::vector< std::vector<int> > two) {
        int count = val.size();

        for (int u = 0; u < count; u++){
            for (int v = 0; v < count; v++) {
                val[u][v] = INT16_MAX - 1;
                for (int one_z = 0; one_z < count; one_z++) {
                    int two_z = one_z;
                    int three_z = one_z;
                    int ones = 0;
                    for (int i = one_z, j = 2, it = 0; it < level; i /= 3, j *= 3, it++) {
                        int mod = i % 3;

                        if (mod == 0) {
                            two_z += j;
                        }

                        if (mod == 1) {
                            ones++;
                        }

                        if (mod == 2) {
                            two_z -= j;
                            three_z -= j;
                        }

                    }

                    val[u][v] = std::min(val[u][v], one[u][one_z] + two[three_z][v] - ones);
                    val[u][v] = std::min(val[u][v], one[u][three_z] + two[one_z][v] - ones);
                    val[u][v] = std::min(val[u][v], one[u][one_z] + two[two_z][v] - ones);
                }
            }
        }
    }

    template<typename Graph>
    void adjust(Graph& g, std::set<std::pair<int, int> >& ae) {
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

    void contract(dominating_set& two, Graph& g, std::set<std::pair<int, int> >& ae) {
        int count = val.size();

        for (int u = 0; u < count; u++) {
            for (int v = 0; v < count; v++) {
                val[u][v] = std::min(two.val[(u * 3) + 1][(v * 3) + 1], two.val[u * 3][v * 3]);

                if ((u % 3 == 1 && check_for_edge(label.first, two.label.first, g, ae)) ||
                (v % 3 == 1 && check_for_edge(label.second, two.label.first, g, ae))) {
                    val[u][v] = std::min(val[u][v], two.val[(u * 3) + 2][(v * 3) + 2]);
                }
            }
        }
    }

    template<typename Graph>
    dominating_set extend(int z, Graph& g, std::set< std::pair<int, int> >& ae) {
        int count = val.size();

        dominating_set res(level + 1);

        std::vector<int> lb;
        std::vector<int> rb;

        get_left_boundary(lb);
        get_right_boundary(rb);

        for (int u = 0; u < count; u++) {
            for (int v = 0; v < count; v++) {
                res.val[u * 3][(v * 3) + 1]= INT16_MAX - 1;
                res.val[(u * 3) + 1][v * 3] = INT16_MAX - 1;
                res.val[(u * 3) + 1][(v * 3) + 2] = INT16_MAX - 1;
                res.val[(u * 3) + 2][(v * 3) + 1] = INT16_MAX - 1;

                res.val[(u * 3) + 2][(v * 3) + 2] = val[u][v];

                if (((u % 3) == 1 && check_for_edge(lb[0], z, g, ae))
                    || ((v % 3) == 1 && check_for_edge(rb[0], z, g, ae))){
                    res.val[u * 3][v * 3] = val[u][v];
                } else {
                    res.val[u * 3][v * 3] = INT16_MAX - 1;
                }

                res.val[u * 3][(v * 3) + 2]= res.val[(u * 3) + 2][(v * 3) + 2];
                res.val[(u * 3) + 2][v * 3] = res.val[(u * 3) + 2][(v * 3) + 2];


                res.val[(u * 3) + 1][(v * 3) + 1] = val[u][v];
                int u_temp = u;
                int v_temp = v;

                if ((u % 3) == 0 && check_for_edge(lb[0], z, g, ae)) {
                    u_temp += 2;
                    res.val[(u * 3) + 1][(v * 3) + 1] = std::min(res.val[(u * 3) + 1][(v * 3) + 1], val[u_temp][v]);
                }

                if ((v % 3) == 0 && check_for_edge(rb[0], z, g, ae)) {
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
    void create(int child_num, Graph& g, std::set<std::pair<int, int> >& ae) {
        const std::vector<int>& children = my_tree->enclosing_tree->t[my_tree->enclosing_face].children;
        std::vector<int> vertices;

        if (child_num < children.size()) {
            dominating_set& child = my_tree->enclosing_tree->t[children[child_num]];
            child.get_left_boundary(vertices);
        } else {
            dominating_set& child = my_tree->enclosing_tree->t[children[child_num - 1]];
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

#endif //TECHNIKA_BAKER_PROBLEMS_HPP
