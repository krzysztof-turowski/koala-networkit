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

    RankPairingHeap() noexcept : firstNode(nullptr), size(0) {}

    ~RankPairingHeap() {
        if (!firstNode) return;
        
        Node* current = firstNode;
        firstNode->next = nullptr;
        
        while (current) {
            Node* next = current->next;
            deleteNode(current);
            current = next;
        }
        
        firstNode = nullptr;
        size = 0;
    }

    void push(const Key& key) override {
        Node* newNode = new Node(key);
        size++;
        
        addNodeToRootList(newNode);
    }

    Key peek() const override {
        if (empty()) throw std::runtime_error("Heap is empty!");
        return firstNode->key;
    }

    Key pop() override {
        if (empty()) throw std::runtime_error("Heap is empty!");
        
        Node* minNode = firstNode;
        Key minKey = minNode->key;
        
        if (size == 1) {
            delete firstNode;
            firstNode = nullptr;
            size = 0;
            return minKey;
        }
        
        addRightSpine(firstNode->left);
        
        std::vector<Node*> consolidatedRoots = consolidateBuckets(firstNode, size);
        
        firstNode = nullptr;
        size--;
        
        push(consolidatedRoots);
        
        delete minNode;
        return minKey;
    }

    bool empty() const override {
        return firstNode == nullptr;
    }

    /**
     * Merges the contents of another heap into this one.
     * After the operation, the other heap will be empty.
     * 
     * @param other The heap to be merged into this one
     */
    void meld(RankPairingHeap& other) {
        if (&other == this) return;
        if (other.empty()) return;
        
        if (empty()) {
            firstNode = other.firstNode;
            size = other.size;
            other.firstNode = nullptr;
            other.size = 0;
            return;
        }
        
        firstNode->next->linkLists(other.firstNode);
        
        if (comp(other.firstNode->key, firstNode->key)) {
            firstNode = other.firstNode;
        }
        
        size += other.size;
        other.firstNode = nullptr;
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
            if (comp(node->key, firstNode->key)) {
                firstNode = node;
            }
            return;
        }
        
        detachFromParent(node, parent);
        
        push(node);
        
        recalculateRank(parent);
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
    Node* firstNode;
    unsigned size;
    Compare comp;

    void detachFromParent(Node* node, Node* parent) {
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

    static Node* linkNodes(Node* x, Node* y, Compare& comparator) {
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

    void addNodeToRootList(Node* node) {
        if (empty()) {
            firstNode = node;
            firstNode->next = firstNode;
            return;
        }
        
        node->next = firstNode->next;
        firstNode->next = node;
        
        if (comp(node->key, firstNode->key)) {
            firstNode = node;
        }
    }

    void push(Node* node) {
        node->parent = nullptr;
        addNodeToRootList(node);
    }

    void push(const std::vector<Node*>& nodes) {
        for (const auto& node : nodes) {
            push(node);
        }
    }

    void addRightSpine(Node* node) {
        while (node) {
            node->parent = nullptr;
            Node* right = node->right;
            node->right = nullptr;
            
            node->next = firstNode->next;
            firstNode->next = node;
            
            node = right;
        }
    }

    int getRank(Node* node) const {
        return node ? node->rank : -1;
    }

    void recalculateRank(Node* node) {
        while (node) {
            int oldRank = node->rank;
            
            if (node->parent == nullptr) {
                node->rank = getRank(node->left) + 1;
                break;
            }
            
            int leftRank = getRank(node->left);
            int rightRank = getRank(node->right);

            int newRank = (std::abs(leftRank - rightRank) > 1) ? std::max(leftRank, rightRank) : std::max(leftRank, rightRank) + 1;
            
            if (newRank >= oldRank) break;
            
            node->rank = newRank;
            node = node->parent;
        }
    }

    void deleteNode(Node* node) {
        if (!node) return;
        if (node->left) deleteNode(node->left);
        if (node->right) deleteNode(node->right);
        delete node;
    }

    std::vector<Node*> consolidateBuckets(Node* root, unsigned nodeCount) {
        std::vector<Node*> result;
        
        if (root->next == root) {
            return result;
        }
        
        constexpr double PHI = (1.0 + std::sqrt(5.0)) / 2.0;
        int maxRank = static_cast<int>(std::log(nodeCount) / std::log(PHI)) + 1;
        std::vector<Node*> buckets(maxRank, nullptr);
        
        auto current = root->next;
        
        buckets[current->rank] = current;
        current = current->next;
        
        while (current != root) {
            auto next = current->next;
            auto rank = current->rank;
            
            if (buckets[rank]) {
                auto merged = linkNodes(current, buckets[rank], comp);
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