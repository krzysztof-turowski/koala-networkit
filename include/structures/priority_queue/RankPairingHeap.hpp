/*
 * RankgPairingHeap.hpp
 *
 * Created on: 21.06.2025
 * Author: Jan Kukowski
 */
#pragma once

#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>
#include <numbers>

#include "structures/PriorityQueue.hpp"

namespace Koala {

/**
 * @ingroup priority-queue
 * 
 * The class implements a priority queue using a rank-pairing heap.
 */
template <class Key, class Compare = std::less<Key>>
class RankPairingHeap : public PriorityQueue<Key, Compare> {
 public:
    struct Node;

    RankPairingHeap() noexcept : first_node(nullptr), size(0) {}

    ~RankPairingHeap() {
        if (!first_node) return;

        Node* current = first_node;
        first_node->next = nullptr;

        while (current) {
            Node* next = current->next;
            delete_node(current);
            current = next;
        }

        first_node = nullptr;
        size = 0;
    }

    void push(const Key& key) override {
        Node* new_node = new Node(key);
        size++;

        add_node_to_root_list(new_node);
    }

    Key& top() override {
        if (empty()) throw std::runtime_error("Heap is empty!");
        return first_node->key;
    }

    void pop() override {
        if (empty()) throw std::runtime_error("Heap is empty!");

        Node* min_node = first_node;

        if (size == 1) {
            delete first_node;
            first_node = nullptr;
            size = 0;
            return;
        }

        add_right_spine(first_node->left);

        std::vector<Node*> consolidated_roots = consolidate_buckets(first_node, size);

        first_node = nullptr;
        size--;

        push(consolidated_roots);

        delete min_node;
    }

    bool empty() const override {
        return first_node == nullptr;
    }

    /**
     * Merges the contents of another heap into this one.
     * After the operation, the other heap will be empty.
     * 
     * @param other The heap to be merged into this one
     */
    void merge(RankPairingHeap& other) {
        if (&other == this) return;
        if (other.empty()) return;

        if (empty()) {
            first_node = other.first_node;
            size = other.size;
            other.first_node = nullptr;
            other.size = 0;
            return;
        }

        first_node->next->linkLists(other.first_node);

        if (comp(other.first_node->key, first_node->key)) {
            first_node = other.first_node;
        }

        size += other.size;
        other.first_node = nullptr;
        other.size = 0;
    }

    /**
     * Decreases a node's key value
     * @param node The node to update
     * @param newKey The new (smaller) key value
     */
    void decreaseKey(Node* node, Key newKey) {
        node->key = newKey;
        Node* parent = node->parent;

        if (!parent) {
            if (comp(node->key, first_node->key)) {
                first_node = node;
            }
            return;
        }

        detach_from_parent(node, parent);

        push(node);

        recalculate_rank(parent);
    }

    struct Node {
        Key key;
        int rank;
        Node* next;
        Node* left;
        Node* right;
        Node* parent;

        Node() : rank(0), next(this), left(nullptr),
                 right(nullptr), parent(nullptr) {}
        explicit Node(Key key) : key(key), rank(0), next(this), left(nullptr),
                               right(nullptr), parent(nullptr) {}

        void getChildren(std::vector<Node*>& result) {
            result.push_back(this);
            if (left) left->getChildren(result);
            if (right) right->getChildren(result);
        }

        void linkLists(Node* other) {
            Node* nextNode = this->next;
            this->next = other->next;
            other->next = nextNode;
        }
    };

 private:
    Node* first_node;
    unsigned size;
    Compare comp;

    void detach_from_parent(Node* node, Node* parent) {
        if (parent->left == node) {
            parent->left = node->right;
        } else if (parent->right == node) {
            parent->right = node->right;
        }

        if (node->right) {
            node->right->parent = parent;
        }

        node->right = nullptr;
        node->parent = nullptr;
    }

    static Node* link_nodes(Node* x, Node* y, Compare& comparator) {
        if (!x) return y;
        if (!y) return x;

        if (comparator(y->key, x->key)) {
            std::swap(x, y);
        }

        y->right = x->left;
        if (y->right) {
            y->right->parent = y;
        }

        x->left = y;
        y->parent = x;
        y->next = nullptr;

        x->rank++;
        return x;
    }

    void add_node_to_root_list(Node* node) {
        if (empty()) {
            first_node = node;
            first_node->next = first_node;
            return;
        }

        node->next = first_node->next;
        first_node->next = node;

        if (comp(node->key, first_node->key)) {
            first_node = node;
        }
    }

    void push(Node* node) {
        node->parent = nullptr;
        add_node_to_root_list(node);
    }

    void push(const std::vector<Node*>& nodes) {
        for (const auto& node : nodes) {
            push(node);
        }
    }

    void add_right_spine(Node* node) {
        while (node) {
            node->parent = nullptr;
            Node* right = node->right;
            node->right = nullptr;

            node->next = first_node->next;
            first_node->next = node;

            node = right;
        }
    }

    int get_rank(Node* node) const {
        return node ? node->rank : -1;
    }

    void recalculate_rank(Node* node) {
        while (node) {
            int oldRank = node->rank;

            if (node->parent == nullptr) {
                node->rank = get_rank(node->left) + 1;
                break;
            }

            int leftRank = get_rank(node->left);
            int rightRank = get_rank(node->right);

            int newRank = (std::abs(leftRank - rightRank) > 1)
                ? std::max(leftRank, rightRank)
                : std::max(leftRank, rightRank) + 1;

            if (newRank >= oldRank) break;

            node->rank = newRank;
            node = node->parent;
        }
    }

    void delete_node(Node* node) {
        if (!node) return;
        if (node->left) delete_node(node->left);
        if (node->right) delete_node(node->right);
        delete node;
    }

    std::vector<Node*> consolidate_buckets(Node* root, unsigned nodeCount) {
        std::vector<Node*> result;

        if (root->next == root) {
            return result;
        }

        int maxRank = static_cast<int>(std::log(nodeCount) / std::log(std::numbers::phi)) + 1;
        std::vector<Node*> buckets(maxRank, nullptr);

        auto current = root->next;

        buckets[current->rank] = current;
        current = current->next;

        while (current != root) {
            auto next = current->next;
            auto rank = current->rank;

            if (buckets[rank]) {
                auto merged = link_nodes(current, buckets[rank], comp);
                buckets[rank] = nullptr;
                result.push_back(merged);
            } else {
                buckets[rank] = current;
            }

            current = next;
        }

        for (auto node : buckets) {
            if (node) result.push_back(node);
        }

        return result;
    }
};

}  // namespace Koala
