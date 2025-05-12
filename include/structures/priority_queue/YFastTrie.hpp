#pragma once

#include "structures/PriorityQueue.hpp"
#include <structures/priority_queue/XFastTrie.hpp>
#include <map>
#include <set>
#include <memory>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <cmath>
#include <random>

namespace Koala {

/**
 * @ingroup priority-queue
 * 
 * The class implements a priority queue using a Y-fast trie tree.
 */
template <class Key, class Compare = std::less<Key>>
class YFastTrie : public PriorityQueue<Key, Compare> {
    static_assert(std::is_unsigned_v<Key>, "YFastTrie only supports unsigned integer keys.");

public:
    explicit YFastTrie(Key universeSize)
        : universeSize(universeSize),
          bitWidth(calculateBitWidth(universeSize)),
          comparator(),
          representatives(universeSize),
          randomEngine(std::random_device{}()),
          distribution(0, 1),
          INVALID_KEY(universeSize + 1),
          minimum(INVALID_KEY) {}

    void push(const Key& key) override {
        if (key >= universeSize) {
            throw std::out_of_range("Key is out of range");
        }
        if (contains(key)) {
            return;
        }
        insertKey(key);
    }

    Key peek() const override {
        if (empty()) {
            throw std::runtime_error("Priority queue is empty");
        }
        return minimum;
    }

    Key pop() override {
        if (empty()) {
            throw std::runtime_error("Priority queue is empty");
        }
        Key minKey = minimum;
        removeKey(minKey);
        return minKey;
    }

    bool empty() const override {
        return size == 0;
    }

private:
    const Key universeSize;
    const int bitWidth;
    Compare comparator;
    XFastTrie<Key, Compare> representatives;
    std::map<Key, std::shared_ptr<std::set<Key>>> buckets;
    size_t size = 0;
    mutable std::mt19937 randomEngine;
    mutable std::uniform_int_distribution<int> distribution;
    const Key INVALID_KEY;
    Key minimum;

    int calculateBitWidth(Key universeSize) {
        return universeSize == 0 ? 0 : static_cast<int>(std::log2(universeSize)) + 1;
    }

    bool contains(Key key) const {
        Key representative = findRepresentative(key);
        if (representative == INVALID_KEY) {
            return false;
        }
        auto it = buckets.find(representative);
        if (it == buckets.end()) {
            return false;
        }
        return it->second->count(key) > 0;
    }

    void insertKey(const Key& key) {
        Key representative = findRepresentative(key);
        std::shared_ptr<std::set<Key>> bucket;
    
        if (representative == INVALID_KEY) {
            bucket = std::make_shared<std::set<Key>>();
            representatives.push(key);
            buckets[key] = bucket;
            representative = key;
        } else {
            bucket = buckets[representative];
        }
    
        bucket->insert(key);
        ++size;

        if (key < minimum) {
            minimum = key;
        }
    
        if (*bucket->begin() != representative) {
            Key newRepresentative = *bucket->begin();
            representatives.remove(representative);
            representatives.push(newRepresentative);
            buckets.erase(representative);
            buckets[newRepresentative] = bucket;
        }
    
        // Randomized bucket split: Split 1/bitWidth of the time
        if (bucket->size() >= 2) {
            double randomValue = distribution(randomEngine);
            if (randomValue < 1.0 / static_cast<double>(bitWidth)) {
                splitBucket(bucket);
            }
        }
    }

    void removeKey(const Key& key) {
        Key representative = findRepresentative(key);
        if (representative == INVALID_KEY) {
            return;
        }

        auto it = buckets.find(representative);
        if (it == buckets.end()) {
            return;
        }

        auto& bucket = it->second;
        if (bucket->erase(key)) {
            --size;

            if (bucket->empty()) {
                representatives.remove(representative);
                buckets.erase(representative);
            } else if (key == representative) {
                Key newRepresentative = *bucket->begin();
                representatives.remove(representative);
                representatives.push(newRepresentative);
                buckets[newRepresentative] = bucket;
                buckets.erase(representative);
            }

            if (key == minimum) {
                minimum = INVALID_KEY;
                for (const auto& [rep, bucket] : buckets) {
                    if (!bucket->empty()) {
                        minimum = *bucket->begin();
                        break;
                    }
                }
            }
        }
    }

    Key findRepresentative(const Key& key) const {
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

        return INVALID_KEY;
    }

    void splitBucket(std::shared_ptr<std::set<Key>>& bucket) {
        auto midIt = bucket->begin();
        std::advance(midIt, bucket->size() / 2);

        std::shared_ptr<std::set<Key>> newBucket = std::make_shared<std::set<Key>>(midIt, bucket->end());
        Key newRepresentative = *newBucket->begin();
        representatives.push(newRepresentative);
        buckets[newRepresentative] = newBucket;

        bucket->erase(midIt, bucket->end());
    }
};

}  // namespace Koala
