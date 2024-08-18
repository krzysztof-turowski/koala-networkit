#pragma once

#include <networkit/Globals.hpp>
#include <vector>

namespace Koala {

/**
 * Implementation of a binary heap with ability to remove elements.
 * Allows for iterating over all elements (not in order).
*/
template<typename Key>
class HeapWithRemove {
 public:
    using handle_type = NetworKit::index;
    using const_iterator = typename std::vector<Key>::const_iterator;
    using iterator = const_iterator;

    HeapWithRemove() {}

    bool empty() const { return heap.empty(); }

    void clear() {
        heap.clear();
        handle.clear();
        heap_index.clear();
    }

    NetworKit::index size() const { return heap.size(); }

    handle_type push(Key value) {
        handle_type ref = heap_index.size();
        NetworKit::index index = heap.size();
        heap.push_back(value);
        handle.push_back(ref);
        heap_index.push_back(index);
        push_up(index);

        return ref;
    }

    Key top() const { return heap.front(); }

    void pop() { erase(handle[0]); }

    void erase(handle_type ref) {
        NetworKit::index index = heap_index[ref];
        NetworKit::index last_index = heap.size() - 1;
        if (index != last_index)
            swap(index, last_index);
        heap.pop_back();
        handle.pop_back();
        if (index != last_index) {
            push_down(index);
            push_up(index);
        }
    }

    void for_elements(const std::function<void(Key)>& handle) const {
        for (auto v : heap)
            handle(v);
    }

    const_iterator begin() const { return heap.begin(); }

    const_iterator end() const { return heap.end(); }

 private:
    std::vector<Key> heap;
    std::vector<handle_type> handle;
    std::vector<NetworKit::index> heap_index;

    NetworKit::index parent_index(NetworKit::index index) { return (index - 1) / 2; }
    NetworKit::index left_son_index(NetworKit::index index) { return 2 * index + 1; }
    NetworKit::index right_son_index(NetworKit::index index) { return 2 * index + 2; }

    void push_up(NetworKit::index index) {
        while (index > 0) {
            NetworKit::index parent = parent_index(index);
            if (heap[index] < heap[parent]) {
                swap(index, parent);
                index = parent;
            } else {
                break;
            }
        }
    }

    void push_down(NetworKit::index index) {
        NetworKit::index minimum = index;
        NetworKit::index left_son = left_son_index(index);
        NetworKit::index right_son = right_son_index(index);;

        if (left_son < heap.size() && heap[left_son] < heap[minimum])
            minimum = left_son;
        if (right_son < heap.size() && heap[right_son] < heap[minimum])
            minimum = right_son;

        if (minimum != index) {
            swap(index, minimum);
            push_down(minimum);
        }
    }

    void swap(NetworKit::index index_1, NetworKit::index index_2) {
        std::swap(heap[index_1], heap[index_2]);
        std::swap(heap_index[handle[index_1]], heap_index[handle[index_2]]);
        std::swap(handle[index_1], handle[index_2]);
    }
};

} /* namespace Koala */
