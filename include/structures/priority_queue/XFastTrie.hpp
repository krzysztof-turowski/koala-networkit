/*
 * XFastTrie.hpp
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
 * The class implements a priority queue using a X-fast trie.
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
        std::shared_ptr<Node> linked_nodes[2] = {nullptr, nullptr};
        std::optional<Key> key;
    };

 public:
    explicit XFastTrie(Key universeSize)
        : universe_size(universeSize),
          bit_width(calculate_bit_width(universeSize)),
          root(std::make_shared<Node>()),
          dummy(std::make_shared<Node>()),
          prefix_levels(bit_width + 1) {
        dummy->linked_nodes[PREV] = dummy;
        dummy->linked_nodes[NEXT] = dummy;
    }

    bool empty() const override {
        return size == 0;
    }

    void push(const Key& key) override {
        if (key >= universe_size) {
            throw std::out_of_range("Key is out of range");
        }
        if (insert(key)) {
            ++size;
        }
    }

    Key& top() override {
        if (empty()) {
            throw std::runtime_error("Priority queue is empty");
        }
        return dummy->linked_nodes[NEXT]->key.value();
    }

    void pop() override {
        if (empty()) {
            throw std::runtime_error("Priority queue is empty");
        }
        Key minKey = dummy->linked_nodes[NEXT]->key.value();
        remove(minKey);
    }

    /**
     * Finds the largest key in the trie that is less than or equal to the given key.
     * 
     * @param key The query key for which to find the predecessor.
     * @return std::optional<Key> The predecessor key if it exists; std::nullopt otherwise.
     */
    std::optional<Key> predecessor(Key key) const {
        auto node = find_predecessor_node(key);
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

        for (int level = 0; level < bit_width; ++level) {
            int bit = get_bit(key, level);
            if (!node->children[bit]) {
                return false;
            }
            node = node->children[bit];
        }

        cleanup_path(node, key);
        unlink_node(node);
        size--;

        return true;
    }

 private:
    Key universe_size;
    int bit_width;
    size_t size = 0;
    std::shared_ptr<Node> root;
    std::shared_ptr<Node> dummy;
    std::vector<std::unordered_map<Key, std::shared_ptr<Node>>> prefix_levels;

    int calculate_bit_width(Key universeSize) {
        return universeSize == 0 ? 1 : static_cast<int>(std::ceil(std::log2(universeSize + 1)));
    }

    int get_bit(Key key, int position) const {
        return (key >> (bit_width - position - 1)) & 1;
    }

    Key get_prefix(Key key, int level) const {
        return key >> (bit_width - level);
    }

    std::shared_ptr<Node> find_predecessor_node(Key key) const {
        int low = 0, high = bit_width + 1;
        auto node = root;

        while (high - low > 1) {
            int mid = (low + high) / 2;
            auto it = prefix_levels[mid].find(get_prefix(key, mid));
            if (it == prefix_levels[mid].end()) {
                high = mid;
            } else {
                node = it->second;
                low = mid;
            }
        }

        if (low == bit_width && node->key.has_value()) {
            return node;
        }

        int dir = get_bit(key, low);
        return (dir == RIGHT)
            ? node->jump
            : (node->jump ? node->jump->linked_nodes[PREV] : nullptr);
    }

    bool insert(Key key) {
        auto node = root;
        int level = 0;

        for (; level < bit_width; ++level) {
            int bit = get_bit(key, level);
            if (!node->children[bit]) break;
            node = node->children[bit];
        }

        if (level == bit_width) {
            return false;
        }

        auto [predecessor, successor] = get_insert_neighbors(node, key, level);

        auto insertedNode = create_path(node, key, level);

        link_neighbors(insertedNode, predecessor, successor);
        update_jump_pointers(insertedNode, key);
        update_prefix_levels(key);

        return true;
    }

    std::pair<std::shared_ptr<Node>, std::shared_ptr<Node>> get_insert_neighbors(
        std::shared_ptr<Node> node,
        Key key,
        int level)
    const {
        int bit = get_bit(key, level);
        if (bit == RIGHT) {
            auto pred = node->jump;
            auto succ = pred ? pred->linked_nodes[NEXT] : dummy->linked_nodes[NEXT];
            return {pred, succ};
        } else {
            auto succ = node->jump;
            auto pred = succ ? succ->linked_nodes[PREV] : dummy->linked_nodes[PREV];
            return {pred, succ};
        }
    }

    std::shared_ptr<Node> create_path(std::shared_ptr<Node> node, Key key, int startLevel) {
        for (int level = startLevel; level < bit_width; ++level) {
            int bit = get_bit(key, level);
            node->children[bit] = std::make_shared<Node>();
            node->children[bit]->parent = node;
            node = node->children[bit];
        }
        node->key = key;
        return node;
    }

    void link_neighbors(
        std::shared_ptr<Node> node,
        std::shared_ptr<Node> pred,
        std::shared_ptr<Node> succ
    ) {
        node->linked_nodes[PREV] = pred;
        node->linked_nodes[NEXT] = succ;
        if (pred) pred->linked_nodes[NEXT] = node;
        if (succ) succ->linked_nodes[PREV] = node;
    }

    void update_jump_pointers(std::shared_ptr<Node> node, Key key) {
        auto parent = node->parent;
        while (parent) {
            if ((parent->children[LEFT] == nullptr && (!parent->jump || parent->jump->key > key)) ||
                (parent->children[RIGHT] == nullptr &&
                (!parent->jump || parent->jump->key < key))) {
                parent->jump = node;
            }
            if ((parent->children[LEFT] != nullptr) &&  (parent->children[RIGHT] != nullptr)) {
                parent->jump = nullptr;
            }
            parent = parent->parent;
        }
    }

    void update_prefix_levels(Key key) {
        auto node = root;
        for (int i = 0; i <= bit_width; ++i) {
            prefix_levels[i][get_prefix(key, i)] = node;
            if (i < bit_width) {
                int bit = get_bit(key, i);
                node = node->children[bit];
            }
        }
    }

    void unlink_node(std::shared_ptr<Node> node) {
        if (node->linked_nodes[PREV]) {
            node->linked_nodes[PREV]->linked_nodes[NEXT] = node->linked_nodes[NEXT];
        }
        if (node->linked_nodes[NEXT]) {
            node->linked_nodes[NEXT]->linked_nodes[PREV] = node->linked_nodes[PREV];
        }
    }

    void cleanup_path(std::shared_ptr<Node> node, Key key) {
        auto parent = node->parent;
        int level = bit_width - 1;

        for (; level >= 0; --level) {
            int bit = get_bit(key, level);

            parent->children[bit] = nullptr;
            prefix_levels[level + 1].erase(get_prefix(key, level + 1));
            if (!parent->children[1 - bit]) parent->jump = nullptr;

            if (parent->children[1 - bit]) break;
            parent = parent->parent;
        }

        if (parent) parent->jump = node;

        for (; level >= 0; --level) {
            if (parent->jump == node) {
                if (!parent->children[LEFT]) {
                    parent->jump = node->linked_nodes[NEXT];
                } else if (!parent->children[RIGHT]) {
                    parent->jump = node->linked_nodes[PREV];
                }
            }
            parent = parent->parent;
        }
    }
};

}  // namespace Koala
