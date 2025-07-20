/*
 * YFastTrie.hpp
 *
 * Created on: 11.05.2025
 * Author: Jan Kukowski
 */

#pragma once

#include "structures/PriorityQueue.hpp"
#include <structures/priority_queue/XFastTrie.hpp>
#include <structures/heap/Treap.hpp>
#include <map>
#include <memory>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <cmath>
#include <optional>

namespace Koala {

/**
 * @ingroup priority-queue
 * 
 * The class implements a priority queue using a Y-fast trie.
 */
template <class Key, class Compare = std::less<Key>>
class YFastTrie : public PriorityQueue<Key, Compare> {
    static_assert(std::is_unsigned_v<Key>, "YFastTrie only supports unsigned integer keys.");

 public:
    explicit YFastTrie(Key universeSize)
        : universe_size(universeSize),
          bit_width(calculate_bit_width(universeSize)),
          representatives(universeSize),
          minimum(std::nullopt) {}

    void push(const Key& key) override {
        if (key >= universe_size) {
            throw std::out_of_range("Key is out of range");
        }
        if (contains(key)) {
            return;
        }
        insert(key);
    }

    Key& top() override {
        if (empty()) {
            throw std::runtime_error("Priority queue is empty");
        }
        return *minimum;
    }

    void pop() override {
        if (empty()) {
            throw std::runtime_error("Priority queue is empty");
        }
        Key minKey = *minimum;
        remove(minKey);
    }

    bool empty() const override {
        return size == 0;
    }

 private:
    Key universe_size;
    int bit_width;
    XFastTrie<Key, Compare> representatives;
    std::optional<Key> minimum;
    size_t size = 0;
    std::map<Key, std::shared_ptr<Koala::Treap<Key>>> buckets;

    int calculate_bit_width(Key universeSize) {
        return universeSize == 0 ? 1 : static_cast<int>(std::log2(universeSize)) + 1;
    }

    bool contains(Key key) const {
        auto representative = find_representative(key);
        if (!representative.has_value()) {
            return false;
        }
        auto it = buckets.find(*representative);
        if (it == buckets.end()) {
            return false;
        }
        return it->second->contains(key);
    }

    void insert(const Key& key) {
        auto representative = find_representative(key);
        std::shared_ptr<Koala::Treap<Key>> bucket;

        if (!representative.has_value()) {
            if (!empty()) {
                Key rep = representatives.top();
                bucket = buckets[rep];
                representative = rep;
            } else {
                bucket = std::make_shared<Koala::Treap<Key>>();
                representatives.push(key);
                buckets[key] = bucket;
                representative = key;
            }
        } else {
            bucket = buckets[*representative];
        }

        bucket->insert(key);
        ++size;

        if (!minimum.has_value() || key < *minimum) {
            minimum = key;
        }

        if (bucket->kth(1) != *representative) {
            Key newRepresentative = bucket->kth(1);
            representatives.remove(*representative);
            representatives.push(newRepresentative);
            buckets.erase(*representative);
            buckets[newRepresentative] = bucket;
            representative = newRepresentative;
        }

        if (bucket->size() > 2 * bit_width) {
            split_bucket(bucket);
        }
    }

    void remove(const Key& key) {
        auto representative = find_representative(key);
        if (!representative.has_value()) {
            return;
        }

        auto it = buckets.find(*representative);

        auto& bucket = it->second;
        if (bucket->contains(key)) {
            bucket->erase(key);
            --size;

            if (bucket->size() == 0) {
                representatives.remove(*representative);
                buckets.erase(*representative);
            } else if (key == *representative) {
                Key newRepresentative = bucket->kth(1);
                representatives.remove(*representative);
                representatives.push(newRepresentative);
                buckets[newRepresentative] = bucket;
                buckets.erase(*representative);
                representative = newRepresentative;
            }

            if (minimum.has_value() && key == *minimum) {
                if (!representatives.empty())
                    minimum = representatives.top();
                else
                    minimum = std::nullopt;
            }
        }

        if (bucket->size() > 0 && bucket->size() < (bit_width + 1) / 2) {
            Key currentRep = bucket->kth(1);
            auto it = buckets.find(currentRep);
            if (it != buckets.end()) {
                auto nextIt = std::next(it);
                if (nextIt != buckets.end()) {
                    auto nextRepresentative = nextIt->first;
                    it->second->merge(*nextIt->second);
                    representatives.remove(nextRepresentative);
                    buckets.erase(nextRepresentative);
                    if (it->second->size() > 2 * bit_width) {
                        split_bucket(it->second);
                    }
                } else if (it != buckets.begin()) {
                    auto prevIt = std::prev(it);
                    prevIt->second->merge(*it->second);
                    representatives.remove(currentRep);
                    buckets.erase(currentRep);
                    if (prevIt->second->size() > 2 * bit_width) {
                        split_bucket(prevIt->second);
                    }
                }
            }
        }
    }

    std::optional<Key> find_representative(const Key& key) const {
        auto it = buckets.find(key);
        if (it != buckets.end()) {
            return key;
        }

        std::optional<Key> predecessor = representatives.predecessor(key);

        if (predecessor.has_value()) {
            auto predecessorIt = buckets.find(predecessor.value());
            if (predecessorIt != buckets.end()) {
                return predecessor.value();
            }
        }

        return std::nullopt;
    }

    void split_bucket(std::shared_ptr<Koala::Treap<Key>>& bucket) {
        auto left = std::make_shared<Koala::Treap<Key>>();
        auto right = std::make_shared<Koala::Treap<Key>>();
        bucket->splitInHalf(*left, *right);

        Key newRepresentative = right->kth(1);
        representatives.push(newRepresentative);
        buckets[newRepresentative] = right;

        Key oldRepresentative = left->kth(1);
        buckets[oldRepresentative] = left;
    }
};

}  // namespace Koala
