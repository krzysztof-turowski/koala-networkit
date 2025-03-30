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

namespace Koala {

/**
 * @ingroup priority-queue
 * 
 * The class implements a priority queue using a Van Emde Boas tree.
 */
class VanEmdeBoasTree : public PriorityQueue<int> {
public:
    VanEmdeBoasTree(int universeSize)
        : universeSize(universeSize), 
          clusterSize(static_cast<int>(std::ceil(std::sqrt(universeSize)))),
          minValue(-1), maxValue(-1) {
        clusters.resize(clusterSize);
    }

    int pop() override {
        if (isEmpty()) {
            throw std::runtime_error("Priority queue is empty");
        }
        int extractedMin = minValue;
        remove(minValue);
        return extractedMin;
    }

    int peek() const override {
        if (isEmpty()) {
            throw std::runtime_error("Priority queue is empty");
        }
        return minValue;
    }

    void push(const int& key) override {
        if (key < 0 || key >= universeSize) {
            throw std::out_of_range("Key is out of range");
        }
        insert(key);
    }

    bool empty() const override {
        return isEmpty();
    }

private:
    int universeSize;
    int clusterSize;
    int minValue, maxValue;
    std::optional<std::unique_ptr<VanEmdeBoasTree>> summaryTree;
    std::vector<std::optional<std::unique_ptr<VanEmdeBoasTree>>> clusters;

    /**
     * Returns the high bits (cluster index) of a key.
     */
    int high(int x) const {
        return x / clusterSize;
    }

    /**
     * Returns the low bits (position within cluster) of a key.
     */
    int low(int x) const {
        return x % clusterSize;
    }

    /**
     * Combines high and low to form a key.
     */
    int index(int high, int low) const {
        return high * clusterSize + low;
    }

    bool isEmpty() const {
        return minValue == -1;
    }

    void insert(int x) {
        if (isEmpty()) {
            minValue = maxValue = x;
            return;
        }
        
        if (x < minValue) {
            std::swap(x, minValue);
        }
        
        if (universeSize > 2) {
            int clusterIndex = high(x);
            int position = low(x);
            
            if (!clusters[clusterIndex].has_value()) {
                clusters[clusterIndex] = std::make_unique<VanEmdeBoasTree>(clusterSize);
            }
            
            if (!summaryTree.has_value()) {
                summaryTree = std::make_unique<VanEmdeBoasTree>(clusterSize);
            }
            
            if (clusters[clusterIndex].value()->isEmpty()) {
                summaryTree.value()->insert(clusterIndex);
                clusters[clusterIndex].value()->minValue = clusters[clusterIndex].value()->maxValue = position;
            } else {
                clusters[clusterIndex].value()->insert(position);
            }
        }
        
        if (x > maxValue) {
            maxValue = x;
        }
    }

    void remove(int x) {
        if (isEmpty()) {
            return;
        }
        
        if (minValue == maxValue) {
            minValue = maxValue = -1;
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
        
        if (x == minValue) {
            int firstCluster = summaryTree.value()->minValue;
            int firstClusterMin = clusters[firstCluster].value()->minValue;

            minValue = index(firstCluster, firstClusterMin);
            clusters[firstCluster].value()->remove(firstClusterMin);
            
            if (clusters[firstCluster].value()->isEmpty()) {
                summaryTree.value()->remove(firstCluster);
            }
        } else {
            int clusterIndex = high(x);
            int position = low(x);
            
            clusters[clusterIndex].value()->remove(position);
            
            if (clusters[clusterIndex].value()->isEmpty()) {
                summaryTree.value()->remove(clusterIndex);
            }
        }
        
        if (x == maxValue) {
            if (summaryTree.value()->isEmpty()) {
                maxValue = minValue;
            } else {
                int lastCluster = summaryTree.value()->maxValue;
                maxValue = index(lastCluster, clusters[lastCluster].value()->maxValue);
            }
        }
    }
};

} // namespace Koala
