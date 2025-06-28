/*
 * VanEmdeBoasTree.hpp
 *
 * Created on: 30.03.2025
 * Author: Jan Kukowski
 */

#pragma once

#include "structures/PriorityQueue.hpp"
#include <cmath>
#include <stdexcept>
#include <vector>
#include <memory>
#include <optional>
#include <type_traits>

namespace Koala {

/**
 * @ingroup priority-queue
 * 
 * The class implements a priority queue using a Van Emde Boas tree.
 */
template <class Key, class Compare = std::less<Key>>
class VanEmdeBoasTree : public PriorityQueue<Key, Compare> {
    static_assert(std::is_unsigned_v<Key>, "VanEmdeBoasTree only supports unsigned integer keys.");

public:
    explicit VanEmdeBoasTree(Key universeSize)
        : universeSize(universeSize), 
          clusterSize(static_cast<Key>(std::ceil(std::sqrt(universeSize)))),
          minValue(std::nullopt), 
          maxValue(std::nullopt) {
        clusters.resize(clusterSize);
    }

    Key pop() override {
        if (isEmpty()) {
            throw std::runtime_error("Priority queue is empty");
        }
        Key extractedMin = *minValue;
        remove(*minValue);
        return extractedMin;
    }

    Key peek() const override {
        if (isEmpty()) {
            throw std::runtime_error("Priority queue is empty");
        }
        return *minValue;
    }

    void push(const Key& key) override {
        if (key >= universeSize) {
            throw std::out_of_range("Key is out of range");
        }
        insert(key);
    }

    bool empty() const override {
        return isEmpty();
    }

private:
    const Key universeSize;
    const Key clusterSize;
    std::optional<Key> minValue, maxValue;
    std::optional<std::unique_ptr<VanEmdeBoasTree>> summaryTree;
    std::vector<std::optional<std::unique_ptr<VanEmdeBoasTree>>> clusters;

    Key high(Key x) const {
        return x / clusterSize;
    }

    Key low(Key x) const {
        return x % clusterSize;
    }

    Key index(Key high, Key low) const {
        return high * clusterSize + low;
    }

    bool isEmpty() const {
        return !minValue.has_value();
    }

    void insert(Key x) {
        if (isEmpty()) {
            minValue = maxValue = x;
            return;
        }

        if (x < *minValue) {
            std::swap(x, *minValue);
        }

        if (universeSize > 2) {
            Key clusterIndex = high(x);
            Key position = low(x);

            if (!clusters[clusterIndex].has_value()) {
                clusters[clusterIndex] = std::make_unique<VanEmdeBoasTree>(clusterSize);
            }

            if (!summaryTree.has_value()) {
                summaryTree = std::make_unique<VanEmdeBoasTree>(clusterSize);
            }

            if (clusters[clusterIndex].value()->isEmpty()) {
                summaryTree.value()->insert(clusterIndex);
                auto &value = clusters[clusterIndex].value();
                value->minValue = value->maxValue = position;
            } else {
                clusters[clusterIndex].value()->insert(position);
            }
        }

        if (x > *maxValue) {
            maxValue = x;
        }
    }

    void remove(Key x) {
        if (isEmpty()) {
            return;
        }

        if (*minValue == *maxValue) {
            minValue = maxValue = std::nullopt;
            return;
        }

        if (universeSize <= 2) {
            if (x == 0) {
                minValue = 1;
            } else {
                minValue = 0;
            }
            maxValue = minValue;
            return;
        }

        if (x == *minValue) {
            Key firstCluster = summaryTree.value()->minValue.value();
            Key firstClusterMin = clusters[firstCluster].value()->minValue.value();

            minValue = index(firstCluster, firstClusterMin);
            clusters[firstCluster].value()->remove(firstClusterMin);

            if (clusters[firstCluster].value()->isEmpty()) {
                summaryTree.value()->remove(firstCluster);
            }
        } else {
            Key clusterIndex = high(x);
            Key position = low(x);

            clusters[clusterIndex].value()->remove(position);

            if (clusters[clusterIndex].value()->isEmpty()) {
                summaryTree.value()->remove(clusterIndex);
            }
        }

        if (x == *maxValue) {
            if (summaryTree.value()->isEmpty()) {
                maxValue = minValue;
            } else {
                Key lastCluster = summaryTree.value()->maxValue.value();
                maxValue = index(lastCluster, clusters[lastCluster].value()->maxValue.value());
            }
        }
    }
};

}  // namespace Koala
