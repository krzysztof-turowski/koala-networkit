/*
 * PriorityQueue.hpp
 *
 * Created on: 29.03.2025
 * Author: Jan Kukowski
 */

#pragma once

#include <functional>

namespace Koala {

/**
 * @ingroup priority-queue
 * The base class for priority queues.
 */
template <class Key, class Compare = std::less<Key>>
class PriorityQueue {
public:
    /**
     * Removes and returns the top element of the priority queue.
     * 
     * @return Key The element with the highest priority.
     * @throws std::runtime_error If the priority queue is empty.
     */
    virtual Key pop() = 0;

    /**
     * Returns the top element of the priority queue without removing it.
     * 
     * @return Key The element with the highest priority.
     * @throws std::runtime_error If the priority queue is empty.
     */
    virtual Key peek() const = 0;

    /**
     * Inserts an element into the priority queue.
     * 
     * @param key The element to insert into the priority queue.
     * @throws std::out_of_range If the key is out of the valid range for the implementation.
     */
    virtual void push(const Key& key) = 0;

    /**
     * Checks if the priority queue is empty.
     * 
     * @return bool True if the priority queue contains no elements, false otherwise.
     */
    virtual bool empty() const = 0;
};

}  // namespace Koala
