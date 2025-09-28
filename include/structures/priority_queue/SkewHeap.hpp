/*
 * SkewHeap.hpp
 *
 * Created on: 19.06.2025
 * Author: Jan Kukowski
 */

#pragma once

#include <functional>
#include <stdexcept>
#include "structures/PriorityQueue.hpp"

namespace Koala {

/**
 * @ingroup priority-queue
 * 
 * The class implements a priority queue using a skew heap.
 */
template <class Key, class Compare = std::less<Key>>
class SkewHeap : public PriorityQueue<Key, Compare> {
 public:
    SkewHeap() noexcept : root(nullptr) {}

    ~SkewHeap() {
        clear(root);
    }

    void push(const Key& key) override {
        Node* new_node = new Node(key);
        root = merge(root, new_node);
    }

    void pop() override {
        if (empty()) {
            throw std::runtime_error("Priority queue is empty");
        }
        Node* old_root = root;
        root = merge(root->left, root->right);
        delete old_root;
    }

    Key& top() override {
        if (empty()) {
            throw std::runtime_error("Priority queue is empty");
        }
        return root->key;
    }

    bool empty() const override {
        return root == nullptr;
    }

 private:
    struct Node {
        Key key;
        Node* left;
        Node* right;
        explicit Node(const Key& k) : key(k), left(nullptr), right(nullptr) {}
    };

    Node* root;
    Compare comp;

    Node* merge(Node* a, Node* b) {
        if (!a) return b;
        if (!b) return a;
        if (comp(b->key, a->key)) {
            std::swap(a, b);
        }
        a->right = merge(a->right, b);
        std::swap(a->left, a->right);
        return a;
    }

    void clear(Node* node) {
        if (!node) return;
        clear(node->left);
        clear(node->right);
        delete node;
    }
};

}  // namespace Koala
