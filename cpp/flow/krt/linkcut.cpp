#include <bits/stdc++.h>

using namespace std;

struct node {
    node *par;
    int id;
    int val;
    set<int> children;
};

struct link_cut_tree {
    vector<node *> _node;

    void init(int n) {
        for (int i = 0; i < n; ++i) {
            node *nn = new node();
            nn->id = i;
            nn->par = nn;
            nn->val = 0;
            _node.push_back(nn);
        }
    }

    void link(int _a, int _b, int v = 0) {
        node *a = _node[_a], *b = _node[_b];
        a->par = b;
        a->val = v;

        b->children.insert(_a);
    }

    void cut(int _a, int _b) {
        node *a = _node[_a], *b = _node[_b];
        a->par = a;

        b->children.erase(_a);
    }

    void modify(int _a, int _b, int v) {
        node *a = _node[_a], *b = _node[_b];

        a->val += v;
        if (a == b) {
            return;
        }
        a = a->par;
        while (a != a->par && a != b) {
            a->val += v;
            a = a->par;
        }
    }

    node *query(int _a, int _b) {
        node *a = _node[_a], *b = _node[_b];
        node *minimum = a;

        if (a == b) {
            return minimum;
        }
        a = a->par;
        while (a != a->par && a != b) {
            if (a->val <= minimum->val) minimum = a;
            a = a->par;
        }
        return minimum;
    }

    node *find_root(int _x) {
        node *x = _node[_x];
        while (x->par != x) {
            x = x->par;
        }
        return x;
    }

    void add_value(int _a, int v) {
        modify(_a, find_root(_a)->id, v);
    }

    node *min(int _x) {
        return query(_x, find_root(_x)->id);
    }

    int get_value(int _x) {
        return query(_x, _x)->val;
    }

    void get_parents() {
        for (int i = 0; i < _node.size(); ++i) {
            if (_node[i]->par != NULL) cout << i << " = " << _node[i]->par->id << endl;
        }
        cout << "END" << endl;
    }

    set<int> getChildren(int _a) {
        return _node[_a]->children;
    }

    node *getParent(int _a) {
        return _node[_a]->par;
    }
};
