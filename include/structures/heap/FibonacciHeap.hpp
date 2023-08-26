/*
 * FibonacciHeap.hpp
 *
 *  Created on: 26.08.2023
 *      Author: Krzysztof Turowski (krzysztof.szymon.turowski@gmail.com)
 */

#pragma once

#include <cassert>
#include <functional>
#include <list>

namespace Koala {

template <class Key, class Compare = std::less<Key> >
class FibonacciHeap {
 private:
		class node {
        friend class FibonacciHeap<Key, Compare>;
        NetworKit::index parent, child, previous, next;
        unsigned flag;
        Key key;

        inline node(NetworKit::index, const Key&);
    };

    std::vector<node> nodes;
    std::list<NetworKit::index> reserved;
    NetworKit::index root;
    std::vector<NetworKit::index> degrees;
		Compare function;

    inline void insert_node(NetworKit::index, NetworKit::index);
    inline void remove_node(NetworKit::index);
 public:
		typedef Key value_type;

		class iterator : public std::iterator<std::input_iterator_tag, NetworKit::index> {
        friend class FibonacciHeap<Key, Compare>;
        NetworKit::index data;
		 public:
        iterator(NetworKit::index data = NetworKit::none) : data(data) { }
        iterator(const iterator &other) : data(other.data) { }
        bool operator==(const iterator &other) { return data == other.data; }
        bool operator!=(const iterator &other) { return !(*this == other); }
        NetworKit::index operator*() const { return data; }
        NetworKit::index* operator->() const { return &**this; }
		};

		inline FibonacciHeap(const Compare& = Compare());

		const Key& top() const;
		iterator push(const Key&);
		void pop();
		void clear();

		void update(iterator, const Key&);
		void erase(iterator);

		void merge(FibonacciHeap&);

		NetworKit::count size() const;
		bool empty() const;

		void check() const;
};

template <class Key, class Compare>
inline FibonacciHeap<Key, Compare>::node::node(NetworKit::index index, const Key &key)
    : parent(NetworKit::none), child(NetworKit::none), previous(index), next(index),
        flag(0), key(key) { }

template <class Key, class Compare>
inline void FibonacciHeap<Key, Compare>::insert_node(NetworKit::index a, NetworKit::index b) {
    auto &A = nodes[a], &B = nodes[b];
    nodes[A.next].previous = B.previous, nodes[B.previous].next = A.next;
    A.next = b, B.previous = a;
}

template <class Key, class Compare>
inline void FibonacciHeap<Key, Compare>::remove_node(NetworKit::index a) {
    auto &A = nodes[a];
    nodes[A.previous].next = A.next, nodes[A.next].previous = A.previous, A.previous = A.next = a;
}

template <class Key, class Compare>
inline FibonacciHeap<Key, Compare>::FibonacciHeap(const Compare &function)
        : root(NetworKit::none), function(function) {
    degrees.resize(sizeof(unsigned) << 3, NetworKit::none);
}

template <class Key, class Compare>
const Key& FibonacciHeap<Key, Compare>::top() const {
    return nodes[root].key;
}

template <class Key, class Compare>
typename FibonacciHeap<Key, Compare>::iterator FibonacciHeap<Key, Compare>::push(const Key &key) {
    NetworKit::index a = nodes.size();
    if (!reserved.empty()) {
        a = reserved.front(), reserved.pop_front();
        nodes[a] = node(a, key);
    } else {
        nodes.push_back(node(a, key));
    }
    if (root == NetworKit::none) {
        root = a;
        return iterator(root);
    }
    insert_node(root, a);
    if (function(nodes[a].key, nodes[root].key)) {
        root = a;
    }
    return iterator(a);
}

template <class Key, class Compare>
void FibonacciHeap<Key, Compare>::pop() {
    NetworKit::index a = nodes[root].child;
    if (a != NetworKit::none) {
        auto b = a;
        do {
            nodes[b].parent = NetworKit::none, b = nodes[b].next;
        } while (a != b);
        insert_node(root, a);
    }
    if (size() == 1) {
        reserved.push_back(root), root = NetworKit::none;
        return;
    }
    unsigned degree_max = 0, degree = 0;
    for (auto a = nodes[root].next, b = nodes[a].next; a != root; degrees[degree] = a, a = b, b = nodes[a].next) {
        while (degrees[degree = nodes[a].flag >> 1] != NetworKit::none) {
            auto c = degrees[degree];
            if (function(nodes[c].key, nodes[a].key)) {
                c = a, a = degrees[degree];
            }
            degrees[degree] = NetworKit::none, remove_node(c), nodes[c].parent = a, nodes[c].flag &= ~1;
            if (nodes[a].child != NetworKit::none) {
                nodes[a].flag += 2, insert_node(nodes[a].child, c);
            } else {
                nodes[a].flag = 2, nodes[a].child = c;
            }
        }
        if (degree > degree_max) {
            degree_max = degree;
        }
    }
    remove_node(root), reserved.push_back(root);

    for (degree = 0; degree <= degree_max; degree++) {
        if (degrees[degree] != NetworKit::none) {
            root = degrees[degree], degrees[degree] = NetworKit::none;
            break;
        }
    }
    for (; degree <= degree_max; degree++) {
        if (degrees[degree] != NetworKit::none) {
            if (function(nodes[degrees[degree]].key, nodes[root].key)) {
                root = degrees[degree];
            }
            degrees[degree] = NetworKit::none;
        }
    }
}

template <class Key, class Compare>
void FibonacciHeap<Key, Compare>::clear() {
    root = NetworKit::none, nodes.clear(), reserved.clear();
}

template <class Key, class Compare>
void FibonacciHeap<Key, Compare>::update(
        typename FibonacciHeap<Key, Compare>::iterator it, const Key &key) {
    NetworKit::index a = *it;
    assert(function(key, nodes[a].key));

    nodes[a].key = key;
    auto b = nodes[a].parent;
    if (b == NetworKit::none) {
        if (function(key, nodes[root].key)) {
            root = a;
        }
        return;
    } else if (!function(key, nodes[b].key)) {
        return;
    }

    while (true) {
        if (nodes[a].next == a) {
            nodes[b].child = NetworKit::none;
        } else {
            if (nodes[b].child == a) {
                nodes[b].child = nodes[a].next;
            }
            remove_node(a), nodes[a].flag &= ~1;
        }
        nodes[b].flag -= 2, insert_node(root, a), nodes[a].parent = NetworKit::none;
        if (function(nodes[a].key, nodes[root].key)) {
            root = a;
        }
        if (nodes[b].parent == NetworKit::none) {
            return;
        }
        if (!(nodes[b].flag & 1)) {
            nodes[b].flag |= 1;
            return;
        }
        a = b, b = nodes[b].parent;
    }
}

template <class Key, class Compare>
void FibonacciHeap<Key, Compare>::erase(iterator it) {
    auto a = *it, b = nodes[a].parent, c = a;
    if (b == NetworKit::none) {
        root = a, pop();
        return;
    }
    while (true) {
        if (a == nodes[a].next) {
            nodes[b].child = NetworKit::none;
        } else {
            if (a == nodes[b].child) {
                nodes[b].child = nodes[a].next;
            }
            remove_node(a), nodes[a].flag &= ~1;
        }
        nodes[b].flag -= 2, insert_node(root, a), nodes[a].parent = NetworKit::none;

        if (nodes[b].parent == NetworKit::none) {
            break;
        }
        if (!(nodes[b].flag & 1)) {
            nodes[b].flag |= 1;
            break;
        }
        a = b, b = nodes[b].parent;
    }
    root = c, pop();
}

template <class Key, class Compare>
void FibonacciHeap<Key, Compare>::merge(FibonacciHeap& heap) {
    if (heap.root == NetworKit::none || this == &heap) {
        return;
    }
    if (root != NetworKit::none) {
        insert_node(root, heap.root);
        if (function(nodes[heap.root].key, nodes[root].key)) {
            root = heap.root;
            // TODO(kturowski): update and merge efficiently reserved
            for (auto &v : heap.reserved) {
                v += size();
            }
            nodes.reserve(nodes.size() + heap.nodes.size());
            nodes.insert(nodes.end(), heap.nodes.begin(), heap.nodes.end());
            reserved.insert(reserved.end(), heap.reserved.begin(), heap.reserved.end());
        }
    } else {
        root = heap.root;
        std::swap(nodes, heap.nodes);
        std::swap(reserved, heap.reserved);
    }
    heap.clear();
}

template <class Key, class Compare>
NetworKit::count FibonacciHeap<Key, Compare>::size() const {
    return nodes.size() - reserved.size();
}

template <class Key, class Compare>
bool FibonacciHeap<Key, Compare>::empty() const {
    return size() == 0;
}

template <class Key, class Compare>
void FibonacciHeap<Key, Compare>::check() const {
    std::queue<NetworKit::index> Q;
    auto a = root;
    do {
        assert(nodes[nodes[a].previous].next == a);
        assert(nodes[nodes[a].next].previous == a);
        assert(nodes[a].parent == NetworKit::none);
        Q.push(a), a = nodes[a].next;
    } while(a != root);

    while (!Q.empty()) {
        a = Q.front(), Q.pop();
        auto child = nodes[a].child;
        NetworKit::count degree = 0;
        if (child != NetworKit::none) {
            do {
                assert(nodes[nodes[child].previous].next == child);
                assert(nodes[nodes[child].next].previous == child);
                assert(nodes[child].parent == a);
                Q.push(child), child = nodes[child].next, degree++;
            } while (child != nodes[a].child);
        }
        assert(degree == (nodes[a].flag >> 1));
    }
}

}  // namespace Koala
