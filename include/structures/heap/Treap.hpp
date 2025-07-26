/*
 * Treap.hpp
 *
 * Created on: 06.06.2025
 * Author: Jan Kukowski
 */

#pragma once

#include <random>
#include <type_traits>
#include <stdexcept>
#include <limits>

namespace Koala {

/**
 * @ingroup heap
 *
 * The class implements a heap using a Treap.
 */
template <class Key, class Compare = std::less<Key>>
class Treap {
    static_assert(std::is_unsigned_v<Key>, "Treap only supports unsigned integer keys.");

 public:
    Treap() noexcept : root(nullptr) {}
    ~Treap() { destroy(root); }

    void insert(Key key) { root = insert(root, key); }
    void erase(Key key) { root = erase(root, key); }
    bool contains(Key key) const { return find(root, key); }
    Key kth(int k) const { return kth(root, k); }
    int size() const { return getSize(root); }

    void merge(Treap& other) {
        root = merge(root, other.root);
        other.root = nullptr;
    }

    void splitInHalf(Treap& leftHalf, Treap& rightHalf) {
        Node *l, *r;
        auto elem = kth(size() / 2);
        split(root, elem, l, r);
        leftHalf.root = l;
        rightHalf.root = r;
        root = nullptr;
    }

 private:
    static std::mt19937 rng;
    static std::uniform_int_distribution<int> dist;

    struct Node {
        Key key;
        int priority;
        int size;
        Node* left;
        Node* right;

        explicit Node(Key k)
            : key(k), priority(dist(rng)), size(1), left(nullptr), right(nullptr) {}
    };

    Node* root;

    static int getSize(Node* t) {
        return t ? t->size : 0;
    }

    static void updateSize(Node* t) {
        if (t) {
            t->size = 1 + getSize(t->left) + getSize(t->right);
        }
    }

    static Node* merge(Node* l, Node* r) {
        if (!l || !r) return l ? l : r;
        if (l->priority > r->priority) {
            l->right = merge(l->right, r);
            updateSize(l);
            return l;
        } else {
            r->left = merge(l, r->left);
            updateSize(r);
            return r;
        }
    }

    static void split(Node* t, Key key, Node*& l, Node*& r) {
        if (!t) {
            l = r = nullptr;
        } else if (t->key <= key) {
            split(t->right, key, t->right, r);
            l = t;
        } else {
            split(t->left, key, l, t->left);
            r = t;
        }
        updateSize(t);
    }

    static bool find(Node* t, Key key) {
        if (!t) return false;
        if (key == t->key) return true;
        return key < t->key ? find(t->left, key) : find(t->right, key);
    }

    static Node* insert(Node* t, Key key) {
        if (find(t, key)) return t;
        Node *l, *r;
        split(t, key, l, r);
        return merge(merge(l, new Node(key)), r);
    }

    static Node* erase(Node* t, Key key) {
        if (!t) return nullptr;
        if (key < t->key) {
            t->left = erase(t->left, key);
        } else if (key > t->key) {
            t->right = erase(t->right, key);
        } else {
            Node* merged = merge(t->left, t->right);
            delete t;
            return merged;
        }
        updateSize(t);
        return t;
    }

    static Key kth(Node* t, int k) {
        if (!t || k <= 0 || k > getSize(t)) {
            throw std::out_of_range("Invalid k");
        }
        int leftSize = getSize(t->left);
        if (k == leftSize + 1) return t->key;
        if (k <= leftSize) return kth(t->left, k);
        return kth(t->right, k - leftSize - 1);
    }

    static void destroy(Node* t) {
        if (!t) return;
        destroy(t->left);
        destroy(t->right);
        delete t;
    }
};

template <class Key, class Compare>
std::mt19937 Treap<Key, Compare>::rng(std::random_device {}());

template <class Key, class Compare>
std::uniform_int_distribution<int> Treap<Key, Compare>::dist(0, std::numeric_limits<int>::max());
}  // namespace Koala
