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
        : universe_size(universeSize),
          cluster_size(static_cast<Key>(std::ceil(std::sqrt(universeSize)))),
          min_value(std::nullopt),
          max_value(std::nullopt) {
        clusters.resize(cluster_size);
    }

    void pop() override {
        if (empty()) {
            throw std::runtime_error("Priority queue is empty");
        }
        remove(*min_value);
    }

    Key& top() override {
        if (empty()) {
            throw std::runtime_error("Priority queue is empty");
        }
        return *min_value;
    }

    void push(const Key& key) override {
        if (key >= universe_size) {
            throw std::out_of_range("Key is out of range");
        }
        insert(key);
    }

    bool empty() const override {
        return !min_value.has_value();
    }

 private:
    const Key universe_size;
    const Key cluster_size;
    std::optional<Key> min_value, max_value;
    std::optional<std::unique_ptr<VanEmdeBoasTree>> summary_tree;
    std::vector<std::optional<std::unique_ptr<VanEmdeBoasTree>>> clusters;

    Key high(Key x) const {
        return x / cluster_size;
    }

    Key low(Key x) const {
        return x % cluster_size;
    }

    Key index(Key high, Key low) const {
        return high * cluster_size + low;
    }

    void insert(Key x) {
        if (empty()) {
            min_value = max_value = x;
            return;
        }

        if (x < *min_value) {
            std::swap(x, *min_value);
        }

        if (universe_size > 2) {
            Key clusterIndex = high(x);
            Key position = low(x);

            if (!clusters[clusterIndex].has_value()) {
                clusters[clusterIndex] = std::make_unique<VanEmdeBoasTree>(cluster_size);
            }

            if (!summary_tree.has_value()) {
                summary_tree = std::make_unique<VanEmdeBoasTree>(cluster_size);
            }

            if (clusters[clusterIndex].value()->empty()) {
                summary_tree.value()->insert(clusterIndex);
                auto &value = clusters[clusterIndex].value();
                value->min_value = value->max_value = position;
            } else {
                clusters[clusterIndex].value()->insert(position);
            }
        }

        if (x > *max_value) {
            max_value = x;
        }
    }

    void remove(Key x) {
        if (empty()) {
            return;
        }

        if (*min_value == *max_value) {
            min_value = max_value = std::nullopt;
            return;
        }

        if (universe_size <= 2) {
            if (x == 0) {
                min_value = 1;
            } else {
                min_value = 0;
            }
            max_value = min_value;
            return;
        }

        if (x == *min_value) {
            Key firstCluster = summary_tree.value()->min_value.value();
            Key firstClusterMin = clusters[firstCluster].value()->min_value.value();

            min_value = index(firstCluster, firstClusterMin);
            clusters[firstCluster].value()->remove(firstClusterMin);

            if (clusters[firstCluster].value()->empty()) {
                summary_tree.value()->remove(firstCluster);
            }
        } else {
            Key clusterIndex = high(x);
            Key position = low(x);

            clusters[clusterIndex].value()->remove(position);

            if (clusters[clusterIndex].value()->empty()) {
                summary_tree.value()->remove(clusterIndex);
            }
        }

        if (x == *max_value) {
            if (summary_tree.value()->empty()) {
                max_value = min_value;
            } else {
                Key lastCluster = summary_tree.value()->max_value.value();
                max_value = index(lastCluster, clusters[lastCluster].value()->max_value.value());
            }
        }
    }
};

}  // namespace Koala
