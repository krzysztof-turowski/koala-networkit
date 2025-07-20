/*
 * WeakHeap.hpp
 *
 * Created on: 14.06.2025
 * Author: Jan Kukowski
 */

#pragma once

#include <vector>
#include <stdexcept>
#include <functional>
#include <cstddef>
#include "structures/PriorityQueue.hpp"

namespace Koala {

/**
 * @ingroup priority-queue
 *
 * The class implements a priority queue using a Weak heap (Dutton, 1993).
 */
template <class Key, class Compare = std::less<Key>>
class WeakHeap : public PriorityQueue<Key, Compare> {
 public:
    WeakHeap() = default;

    void push(const Key& key) override {
        std::size_t index = data.size();
        data.push_back(key);
        flip.push_back(0);

        if ((index % 2 == 0) && index > 0) {
            flip[index / 2] = 0;
        }

        sift_up(index);
    }

    void pop() override {
        if (empty()) {
            throw std::runtime_error("Priority queue is empty");
        }
        std::size_t last = data.size() - 1;
        data[0] = data[last];
        data.pop_back();
        flip.pop_back();
        if (last > 1) {
            sift_down(0);
        }
    }

    Key& top() override {
        if (empty()) {
            throw std::runtime_error("Priority queue is empty");
        }
        return data[0];
    }

    bool empty() const override {
        return data.empty();
    }

 private:
    std::vector<Key> data;
    std::vector<bool> flip;
    Compare comp;

    std::size_t distinguished_ancestor(std::size_t index) const {
        while ((index % 2) == flip[index / 2]) {
            index = index / 2;
        }
        return index / 2;
    }

    bool join(std::size_t parent, std::size_t child) {
        if (comp(data[child], data[parent])) {
            std::swap(data[parent], data[child]);
            flip[child] = !flip[child];
            return false;
        }
        return true;
    }

    void sift_up(std::size_t start) {
        std::size_t current = start;
        while (current != 0) {
            std::size_t ancestor = distinguished_ancestor(current);
            if (join(ancestor, current)) {
                break;
            }
            current = ancestor;
        }
    }

    void sift_down(std::size_t start) {
        std::size_t size = data.size();
        std::size_t descendant = 2 * start + 1 - flip[start];
        while (2 * descendant + flip[descendant] < size) {
            descendant = 2 * descendant + flip[descendant];
        }
        while (descendant != start) {
            join(start, descendant);
            descendant = descendant / 2;
        }
    }
};

}  // namespace Koala
