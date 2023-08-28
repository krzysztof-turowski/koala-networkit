/*
 * BinomialHeap.hpp
 *
 *  Created on: 24.08.2023
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

#include <iterator>
#include <list>

#include <networkit/graph/Graph.hpp>

namespace Koala {

/**
 * @ingroup heap
 * Heap structure from Vuillemin "A Data Structure for Manipulating Priority Queues" (1978).
 */
template <class Key, class Compare = std::less<Key>>
class BinomialHeap {
 private:
    class node {
     public:
        NetworKit::index parent, child, next;
        NetworKit::count degree;
        Key key;

        explicit inline node(const Key &key)
            : parent(NetworKit::none), child(NetworKit::none), next(NetworKit::none),
                degree(0), key(key) { }
    };

    std::vector<node> nodes;
    std::list<NetworKit::index> reserved;
    NetworKit::index root, minimum;
    Compare function;

    void insert_node(NetworKit::index, NetworKit::index);
    inline NetworKit::index join(NetworKit::index, NetworKit::index);
    inline NetworKit::index reverse(NetworKit::index);

 public:
    typedef Key value_type;

    class iterator : public std::iterator<std::input_iterator_tag, NetworKit::index> {
        NetworKit::index data;
     public:
        explicit iterator(NetworKit::index data = NetworKit::none) : data(data) { }
        iterator(const iterator &other) : data(other.data) { }
        bool operator==(const iterator &other) { return data == other.data; }
        bool operator!=(const iterator &other) { return !(*this == other); }
        NetworKit::index operator*() const { return data; }
        NetworKit::index* operator->() const { return &**this; }
    };

    explicit BinomialHeap(const Compare& = Compare());

    const value_type& top() const;
    iterator push(const value_type&);
    void pop();
    void clear();

    void update(iterator, const value_type&);
    void erase(iterator);

    NetworKit::count size() const;
    bool empty() const;

    void check() const;
};

template <class Key, class Compare>
inline void BinomialHeap<Key, Compare>::insert_node(NetworKit::index a, NetworKit::index b) {
    auto &A = nodes[a], &B = nodes[b];
    B.parent = a, B.next = A.child, A.child = b, ++A.degree;
}

template <class Key, class Compare>
inline NetworKit::index BinomialHeap<Key, Compare>::join(NetworKit::index a, NetworKit::index b) {
    if (nodes[a].degree < nodes[b].degree) {
        std::swap(a, b);
    }
    auto start = a, c = a;
    for (a = nodes[a].next; a != NetworKit::none && b != NetworKit::none; a = nodes[a].next) {
        if (nodes[a].degree < nodes[b].degree) {
            std::swap(a, b);
        }
        c = nodes[c].next = a;
    }
    nodes[c].next = a != NetworKit::none ? a : b;

    start = reverse(start);
    a = NetworKit::none, b = start, c = nodes[b].next;
    while (c != NetworKit::none) {
        assert(a == NetworKit::none || nodes[a].next == b);
        if (nodes[b].degree != nodes[c].degree || (nodes[c].next != NetworKit::none
                && nodes[c].degree == nodes[nodes[c].next].degree)) {
            a = b, b = c;
        } else if (!function(nodes[b].key, nodes[c].key)) {
            nodes[b].next = nodes[c].next, insert_node(b, c);
        } else {
            if (a != NetworKit::none) {
                nodes[a].next = c;
            } else {
                start = c;
            }
            insert_node(c, b), b = c;
        }
        c = nodes[b].next;
    }
    return reverse(start);
}

template <class Key, class Compare>
inline NetworKit::index BinomialHeap<Key, Compare>::reverse(NetworKit::index a) {
    auto b = nodes[a].next;
    nodes[a].next = NetworKit::none;
    while (b != NetworKit::none) {
        auto c = nodes[b].next;
        nodes[b].next = a, a = b, b = c;
    }
    return a;
}

template <class Key, class Compare>
BinomialHeap<Key, Compare>::BinomialHeap(const Compare &function)
    : root(NetworKit::none), minimum(NetworKit::none), function(function) { }

template <class Key, class Compare>
const typename BinomialHeap<Key, Compare>::value_type& BinomialHeap<Key, Compare>::top() const {
    return nodes[minimum].key;
}

template <class Key, class Compare>
typename BinomialHeap<Key, Compare>::iterator BinomialHeap<Key, Compare>::push(
        const typename BinomialHeap<Key, Compare>::value_type &key) {
    NetworKit::index a = nodes.size();
    if (!reserved.empty()) {
        a = reserved.front(), reserved.pop_front();
        nodes[a] = node(key);
    } else {
        nodes.push_back(node(key));
    }
    if (root == NetworKit::none) {
        root = minimum = a;
        return iterator(root);
    }
    assert(root != NetworKit::none);
    assert(a != NetworKit::none);
    root = join(root, a);
    if (!function(nodes[a].key, nodes[minimum].key)) {
        minimum = a;
    }
    return iterator(a);
}

template <class Key, class Compare>
void BinomialHeap<Key, Compare>::pop() {
    if (root == minimum) {
        root = nodes[root].next;
    } else {
        auto a = root;
        while (nodes[a].next != minimum) {
            a = nodes[a].next;
        }
        nodes[a].next = nodes[minimum].next;
    }
    if (size() == 1) {
        reserved.push_back(minimum), minimum = NetworKit::none;
        return;
    }
    auto child = nodes[minimum].child;
    if (child != NetworKit::none) {
        for (auto a = child; a != NetworKit::none; a = nodes[a].next) {
            nodes[a].parent = NetworKit::none;
        }
        assert(child != NetworKit::none);
        root = root == NetworKit::none ? child : join(root, child);
    }
    reserved.push_back(minimum), minimum = root;
    if (minimum != NetworKit::none) {
        for (auto a = nodes[root].next; a != NetworKit::none; a = nodes[a].next) {
            if (!function(nodes[a].key, nodes[minimum].key)) {
                minimum = a;
            }
        }
    }
}

template <class Key, class Compare>
void BinomialHeap<Key, Compare>::clear() {
    root = minimum = NetworKit::none, nodes.clear(), reserved.clear();
}

template <class Key, class Compare>
void BinomialHeap<Key, Compare>::update(
        typename BinomialHeap<Key, Compare>::iterator it,
        const typename BinomialHeap<Key, Compare>::value_type &key) {
    auto a = *it;
    assert(!function(key, nodes[a].key));
    nodes[a].key = key;
    if (!function(key, nodes[minimum].key)) {
        minimum = a;
    }

    if (nodes[a].parent == NetworKit::none
            || !function(nodes[nodes[a].parent].key, nodes[a].key)) {
        return;
    }

    auto start = NetworKit::none, previous = NetworKit::none, b = a, c = nodes[a].parent;
    while (c != NetworKit::none) {
        auto d = nodes[c].child;
        nodes[c].child = nodes[b].next;
        if (b == d) {
            d = c, --nodes[c].degree;
        } else {
            auto e = d;
            while (b != nodes[e].next) {
                e = nodes[e].next, --nodes[c].degree;
            }
            c = nodes[e].next, nodes[c].degree -= 2;
        }
        nodes[b].next = start, b = c, c = nodes[c].parent, start = previous, previous = d;
    }
    if (b != root) {
        root = nodes[root].next;
    } else {
        c = root;
        while (b != nodes[c].next) {
            c = nodes[c].next;
        }
        nodes[c].next = nodes[b].next;
    }
    nodes[b].next = start, start = previous;
    if (start != NetworKit::none) {
        for (b = start; b != NetworKit::none; b = nodes[b].next) {
            nodes[b].parent = NetworKit::none;
        }
        root = root ? join(root, start) : start;
    }
    nodes[a].parent = nodes[a].next = NetworKit::none;
    root = root != NetworKit::none ? join(root, a) : a;
}

template <class Key, class Compare>
void BinomialHeap<Key, Compare>::erase(typename BinomialHeap<Key, Compare>::iterator it) {
    auto a = *it;
    if (a == root) {
        pop();
        return;
    }
    auto start = nodes[a].child, previous = nodes[a].child, b = a, c = nodes[a].parent;
    while (c != NetworKit::none) {
        auto next = nodes[c].child;
        nodes[c].child = nodes[b].next;
        if (b == next) {
            next = c, --nodes[c].degree;
        } else {
          auto d = next;
          while (b != nodes[d].next) {
              d = nodes[d].next, --nodes[c].degree;
          }
          nodes[d].next = c, nodes[c].degree -= 2;
        }
        nodes[b].next = start, b = c, c = nodes[c].parent, start = previous, previous = next;
    }
    if (b == root) {
        root = nodes[root].next;
    } else {
        c = root;
        while (b != nodes[c].next) {
            c = nodes[c].next;
        }
        nodes[c].next = nodes[b].next;
    }
    nodes[b].next = start, start = previous;
    if (start != NetworKit::none) {
        for (b = start; b != NetworKit::none; b = nodes[b].next) {
            nodes[b].parent = NetworKit::none;
        }
    }
    root = root != NetworKit::none ? join(root, start) : start;
    if (minimum == a) {
        for (minimum = root, b = nodes[root].next; b != NetworKit::none; b = nodes[b].next) {
            if (!function(nodes[b].key, nodes[minimum].key)) {
                minimum = b;
            }
        }
    }
    reserved.push_back(a);
}

template <class Key, class Compare>
NetworKit::count BinomialHeap<Key, Compare>::size() const {
    return nodes.size() - reserved.size();
}

template <class Key, class Compare>
bool BinomialHeap<Key, Compare>::empty() const {
    return root == NetworKit::none;
}

template <class Key, class Compare>
void BinomialHeap<Key, Compare>::check() const {
    if (empty()) {
        return;
    }
    assert(nodes[minimum].parent == NetworKit::none);
    std::queue<NetworKit::index> Q;
    NetworKit::count total = 0;
    for (auto a = root; a != NetworKit::none; a = nodes[a].next) {
        assert(nodes[a].parent == NetworKit::none);
        Q.push(a);
    }
    while (!Q.empty()) {
        auto a = Q.front();
        Q.pop(), ++total;
        if (nodes[a].next != NetworKit::none) {
            assert(nodes[a].degree > nodes[nodes[a].next].degree);
        }
        auto child = nodes[a].child;
        NetworKit::count degree = 0;
        while (child != NetworKit::none) {
            assert(nodes[child].parent == a);
            assert(!function(nodes[a].key, nodes[child].key));
            Q.push(child), child = nodes[child].next, ++degree;
        }
        assert(degree == nodes[a].degree);
    }
    assert(total == size());
}

}  // namespace Koala
