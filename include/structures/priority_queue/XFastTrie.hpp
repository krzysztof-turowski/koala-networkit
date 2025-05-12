/*
 * VanEmdeBoasTree.hpp
 *
 * Created on: 11.05.2025
 * Author: Jan Kukowski
 */

#pragma once

#include "structures/PriorityQueue.hpp"
#include <memory>
#include <cmath>
#include <limits>
#include <unordered_map>
#include <vector>
#include <optional>

namespace Koala {

/**
 * @ingroup priority-queue
 * 
 * The class implements a priority queue using a X-fast trie tree.
 */
template <class Key, class Compare = std::less<Key>>
class XFastTrie : public PriorityQueue<Key, Compare> {
    static_assert(std::is_unsigned_v<Key>, "XFastTrie only supports unsigned integer keys.");

    static constexpr int LEFT = 0;
    static constexpr int RIGHT = 1;

    static constexpr int PREV = 0;
    static constexpr int NEXT = 1;

    struct Node {
        std::shared_ptr<Node> children[2] = {nullptr, nullptr};
        std::shared_ptr<Node> jump = nullptr;
        std::shared_ptr<Node> parent = nullptr;
        std::shared_ptr<Node> linkedNodes[2] = {nullptr, nullptr};
        std::optional<Key> key;
    };

public:
    explicit XFastTrie(Key universeSize)
        : universeSize(universeSize),
          bitWidth(calculateBitWidth(universeSize)),
          root(std::make_shared<Node>()),
          dummy(std::make_shared<Node>()),
          prefixLevels(bitWidth + 1)
    {
        dummy->linkedNodes[PREV] = dummy;
        dummy->linkedNodes[NEXT] = dummy;
    }

    bool empty() const override {
        return size == 0;
    }

    void push(const Key& key) override {
        if (key >= universeSize) {
            throw std::out_of_range("Key is out of range");
        }
        if (insert(key)) {
            ++size;
        }
    }

    Key peek() const override {
        if (empty()) {
            throw std::runtime_error("Priority queue is empty");
        }
        return dummy->linkedNodes[NEXT]->key.value();
    }

    Key pop() override {
        if (empty()) {
            throw std::runtime_error("Priority queue is empty");
        }
        Key minKey = dummy->linkedNodes[NEXT]->key.value();
        if (remove(minKey)) {
            --size;
        }
        return minKey;
    }
    
    /**
     * Finds the largest key in the trie that is less than or equal to the given key.
     * 
     * @param key The query key for which to find the predecessor.
     * @return std::optional<Key> The predecessor key if it exists; std::nullopt otherwise.
     */
    std::optional<Key> predecessor(Key key) const {
        auto node = findPredecessorNode(key);
        if (!node) {
            return std::nullopt;
        }
        return node->key;
    }

    /**
     * Removes a key from the XFastTrie if it exists.
     * 
     * @param key The key to be removed from the trie.
     * @return true If the key was successfully found and removed.
     * @return false If the key does not exist in the trie.
     */
    bool remove(Key key) {
        auto node = root;

        for (int level = 0; level < bitWidth; ++level) {
            int bit = getBit(key, level);
            if (!node->children[bit]) return false;
            node = node->children[bit];
        }

        unlinkNode(node);
        cleanupPath(node, key);
        erasePrefixes(key);

        return true;
    }

private:
    Key universeSize;
    int bitWidth;
    size_t size = 0;
    std::shared_ptr<Node> root;
    std::shared_ptr<Node> dummy;
    std::vector<std::unordered_map<Key, std::shared_ptr<Node>>> prefixLevels;

    int calculateBitWidth(Key universeSize) {
        return universeSize == 0 ? 1 : static_cast<int>(std::ceil(std::log2(universeSize + 1)));
    }

    int getBit(Key key, int position) const {
        return (key >> (bitWidth - position - 1)) & 1;
    }

    Key getPrefix(Key key, int level) const {
        return key >> (bitWidth - level) >> 1;
    }

    std::shared_ptr<Node> findPredecessorNode(Key key) const {
        int low = 0, high = bitWidth + 1;
        auto node = root;

        while (high - low > 1) {
            int mid = (low + high) / 2;
            auto it = prefixLevels[mid].find(getPrefix(key, mid));
            if (it == prefixLevels[mid].end()) {
                high = mid;
            } else {
                node = it->second;
                low = mid;
            }
        }

        if (low == bitWidth && node->key.has_value()) {
            return node;
        }

        int dir = getBit(key, low);
        return (dir == RIGHT) ? node->jump : (node->jump ? node->jump->children[LEFT] : nullptr);
    }

    bool insert(Key key) {
        auto node = root;
        int level = 0;

        for (; level < bitWidth; ++level) {
            int bit = getBit(key, level);
            if (!node->children[bit]) break;
            node = node->children[bit];
        }

        if (level == bitWidth) {
            return false;
        }

        auto [predecessor, successor] = getInsertNeighbors(node, key, level);

        node->jump = nullptr;
        auto insertedNode = createPath(node, key, level);

        linkNeighbors(insertedNode, predecessor, successor);
        updateJumpPointers(insertedNode, key);
        updatePrefixLevels(key);

        return true;
    }

    std::pair<std::shared_ptr<Node>, std::shared_ptr<Node>> getInsertNeighbors(std::shared_ptr<Node> node, Key key, int level) const {
        int bit = getBit(key, level);
        if (bit == RIGHT) {
            auto pred = node->jump;
            auto succ = pred ? pred->linkedNodes[NEXT] : dummy->linkedNodes[NEXT];
            return {pred, succ};
        } else {
            auto succ = node->jump;
            auto pred = succ ? succ->linkedNodes[PREV] : dummy->linkedNodes[PREV];
            return {pred, succ};
        }
    }

    std::shared_ptr<Node> createPath(std::shared_ptr<Node> node, Key key, int startLevel) {
        for (int level = startLevel; level < bitWidth; ++level) {
            int bit = getBit(key, level);
            node->children[bit] = std::make_shared<Node>();
            node->children[bit]->parent = node;
            node = node->children[bit];
        }
        node->key = key;
        return node;
    }

    void linkNeighbors(std::shared_ptr<Node> node, std::shared_ptr<Node> pred, std::shared_ptr<Node> succ) {
        node->linkedNodes[PREV] = pred;
        node->linkedNodes[NEXT] = succ;
        if (pred) pred->linkedNodes[NEXT] = node;
        if (succ) succ->linkedNodes[PREV] = node;
    }

    void updateJumpPointers(std::shared_ptr<Node> node, Key key) {
        auto parent = node->parent;
        while (parent) {
            if ((parent->children[LEFT] == nullptr && (!parent->jump || parent->jump->key > key)) ||
                (parent->children[RIGHT] == nullptr && (!parent->jump || parent->jump->key < key))) {
                parent->jump = node;
            }
            parent = parent->parent;
        }
    }

    void updatePrefixLevels(Key key) {
        auto node = root;
        for (int i = 0; i <= bitWidth; ++i) {
            prefixLevels[i][getPrefix(key, i)] = node;
            if (i < bitWidth) {
                int bit = getBit(key, i);
                node = node->children[bit];
            }
        }
    }

    void unlinkNode(std::shared_ptr<Node> node) {
        if (node->linkedNodes[PREV]) {
            node->linkedNodes[PREV]->linkedNodes[NEXT] = node->linkedNodes[NEXT];
        }
        if (node->linkedNodes[NEXT]) {
            node->linkedNodes[NEXT]->linkedNodes[PREV] = node->linkedNodes[PREV];
        }
    }

    void cleanupPath(std::shared_ptr<Node> node, Key key) {
        auto parent = node->parent;
        for (int level = bitWidth - 1; level >= 0; --level) {
            int bit = getBit(key, level);
            parent->children[bit] = nullptr;

            if (parent->children[1 - bit]) break;
            parent = parent->parent;
        }

        while (parent) {
            if (parent->jump == node) {
                int bit = getBit(key, 0);
                parent->jump = parent->children[1 - bit] ? parent->children[1 - bit] : nullptr;
            }
            parent = parent->parent;
        }
    }

    void erasePrefixes(Key key) {
        auto node = root;
        for (int i = 0; i <= bitWidth; ++i) {
            prefixLevels[i].erase(getPrefix(key, i));
            if (i < bitWidth) {
                int bit = getBit(key, i);
                if (!node->children[bit]) break;
                node = node->children[bit];
            }
        }
    }
};

} // namespace Koala
