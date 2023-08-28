/*
 * PairingHeap.hpp
 *
 *  Created on: 27.08.2023
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

#include <iterator>

namespace Koala {

/**
 * @ingroup heap
 * Heap structure from Fredman et al. "The Pairing Heap: A New Form of Self-Adjusting Heap" (1986).
 */
template <class Key, class Compare = std::less<Key>>
class PairingHeap {
 private:
    class node {
     public:
        NetworKit::index parent, child, previous, next;
        NetworKit::count degree;
        Key key;

        explicit inline node(const Key &key)
            : parent(NetworKit::none), child(NetworKit::none), previous(NetworKit::none),
                next(NetworKit::none), degree(0), key(key) { }
    };

    std::vector<node> nodes;
    std::list<NetworKit::index> reserved;
    NetworKit::index root;
    Compare function;

    inline void insert_node(NetworKit::index, NetworKit::index);
    inline void remove_node(NetworKit::index);

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

    explicit PairingHeap(const Compare& = Compare());

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
inline void PairingHeap<Key, Compare>::insert_node(NetworKit::index a, NetworKit::index b) {
    auto &A = nodes[a], &B = nodes[b];
    if (A.child != NetworKit::none) {
        nodes[A.child].previous = b;
    }
    B.parent = a, B.previous = NetworKit::none, B.next = A.child, A.child = b, ++A.degree;
}

template <class Key, class Compare>
inline void PairingHeap<Key, Compare>::remove_node(NetworKit::index a) {
    auto &A = nodes[a];
    if (nodes[A.parent].child == a) {
        nodes[A.parent].child = A.next;
    } else {
        nodes[A.previous].next = A.next;
    }
    if (A.next != NetworKit::none) {
        nodes[A.next].previous = A.previous;
    }
    --nodes[A.parent].degree, A.parent = A.previous = A.next = NetworKit::none;
}

template <class Key, class Compare>
inline PairingHeap<Key, Compare>::PairingHeap(const Compare &function)
    : root(NetworKit::none), function(function) { }

template <class Key, class Compare>
const typename PairingHeap<Key, Compare>::value_type& PairingHeap<Key, Compare>::top() const {
    return nodes[root].key;
}

template <class Key, class Compare>
typename PairingHeap<Key, Compare>::iterator PairingHeap<Key, Compare>::push(
        const typename PairingHeap<Key, Compare>::value_type &key) {
    NetworKit::index a = nodes.size();
    if (!reserved.empty()) {
        a = reserved.front(), reserved.pop_front();
        nodes[a] = node(key);
    } else {
        nodes.push_back(node(key));
    }
    if (root == NetworKit::none) {
        root = a;
        return iterator(root);
    }
    if (!function(nodes[a].key, nodes[root].key)) {
        insert_node(a, root), root = a;
    } else {
        insert_node(root, a);
    }
    return iterator(a);
}

template <class Key, class Compare>
void PairingHeap<Key, Compare>::pop() {
    if (size() == 1) {
        reserved.push_back(root), root = NetworKit::none;
        return;
    }

    auto a = nodes[root].child, b = a;
    reserved.push_back(root), root = a;
    while (a != NetworKit::none) {
        b = nodes[a].next;
        if (b == NetworKit::none) {
            nodes[a].parent = NetworKit::none, b = a;
            break;
        }
        if (!function(nodes[a].key, nodes[b].key)) {
            if (nodes[b].next != NetworKit::none) {
                nodes[nodes[b].next].previous = a;
            }
            nodes[a].next = nodes[b].next, nodes[a].parent = NetworKit::none;
            insert_node(a, b), b = a, a = nodes[a].next;
        } else {
            if (nodes[a].previous != NetworKit::none) {
                nodes[nodes[a].previous].next = b;
            }
            nodes[b].previous = nodes[a].previous, nodes[b].parent = NetworKit::none;
            insert_node(b, a), a = nodes[b].next;
        }
    }

    root = b, a = nodes[b].previous;
    while (a != NetworKit::none) {
        if (!function(nodes[a].key, nodes[root].key)) {
            insert_node(a, root), nodes[a].next = NetworKit::none, root = a;
        } else {
            nodes[root].previous = nodes[a].previous, insert_node(root, a);
        }
        a = nodes[root].previous;
    }
}

template <class Key, class Compare>
void PairingHeap<Key, Compare>::clear() {
    root = NetworKit::none, nodes.clear(), reserved.clear();
}

template <class Key, class Compare>
void PairingHeap<Key, Compare>::update(
        typename PairingHeap<Key, Compare>::iterator it,
        const typename PairingHeap<Key, Compare>::value_type &key) {
    NetworKit::index a = *it;
    assert(!function(key, nodes[a].key));
    nodes[a].key = key;
    if (nodes[a].parent == NetworKit::none) {
        return;
    }
    remove_node(a);
    if (!function(nodes[a].key, nodes[root].key)) {
        insert_node(a, root), root = a;
    } else {
        insert_node(root, a);
    }
}

template <class Key, class Compare>
void PairingHeap<Key, Compare>::erase(iterator it) {
    NetworKit::index a = *it;
    if (nodes[a].parent != NetworKit::none) {
        remove_node(a), insert_node(a, root), root = a;
    }
    pop();
}

template <class Key, class Compare>
NetworKit::count PairingHeap<Key, Compare>::size() const {
    return nodes.size() - reserved.size();
}

template <class Key, class Compare>
bool PairingHeap<Key, Compare>::empty() const {
    return root == NetworKit::none;
}

template <class Key, class Compare>
void PairingHeap<Key, Compare>::check() const {
    if (empty()) {
        return;
    }
    std::queue<NetworKit::index> Q;
    Q.push(root);
    assert(nodes[root].previous == NetworKit::none);
    assert(nodes[root].next == NetworKit::none);
    assert(nodes[root].parent == NetworKit::none);
    NetworKit::count total = 0;
    while (!Q.empty()) {
        auto a = Q.front();
        Q.pop(), ++total;
        assert(nodes[a].previous == NetworKit::none || nodes[nodes[a].previous].next == a);
        assert(nodes[a].next == NetworKit::none || nodes[nodes[a].next].previous == a);
        auto child = nodes[a].child;
        assert(child == NetworKit::none || nodes[child].previous == NetworKit::none);
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
