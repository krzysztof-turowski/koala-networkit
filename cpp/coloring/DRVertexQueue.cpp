/*
* DRVertexQueue.cpp
*
*  Created on: 28.02.2023
*      Author: Zofia Glapa (zofia.glapa@student.uj.edu.pl)
*/

#include <coloring/DRVertexQueue.hpp>

namespace Koala {

    bool operator>(const NodeValuePair& a, const NodeValuePair& b) {
        if (a.value == b.value) {
            return a.node < b.node;
        }
        return a.value > b.value;
    }

    bool operator==(const NodeValuePair& a, const NodeValuePair& b) {
        return a.node == b.node && a.value == b.value;
    }

    DRVertexQueue::DRVertexQueue() {
        size = 0;
    }

    DRVertexQueue::DRVertexQueue(const std::vector<std::pair<int, int>> & nodes) {
        size = nodes.size();
        queue.resize(size);
        for (int i = 0; i < size; i++) {
            queue[i] = {nodes[i].first, nodes[i].second};
        }
        for (int i = 0; i < size; i++) {
            position[queue[i].node] = i;
        }
        buildHeap();
    }

    void DRVertexQueue::updateValue(const int & v, const int& value) {
        int i = position[v];
        queue[i].value = value;
        heapify(i/2);
    }

    void DRVertexQueue::insert(const int & v, const int& value) {
        queue.push_back({v, value});
        position[v] = size;
        size++;
        heapify(size - 1);
    }

    NodeValuePair DRVertexQueue::pop() {
        auto v = queue[0];
        swap(0, size - 1);
        size--;
        heapify(0);
        return v;
    }

    bool DRVertexQueue::empty() {
        return size == 0;
    }

    void DRVertexQueue::swap(const int& i, const int& j) {
        std::swap(queue[i], queue[j]);
        std::swap(position[queue[i].node], position[queue[j].node]);
    }

    void DRVertexQueue::heapify(const int& i) {
        int l = 2 * i + 1;
        int r = 2 * i + 2;
        int biggest = i;
        if (l < size && queue[l] > queue[biggest]) {
            biggest = l;
        }
        if (r < size && queue[r] > queue[biggest]) {
            biggest = r;
        }
        if (biggest != i) {
            swap(i, biggest);
            heapify(i/2);
        }
    }

    void DRVertexQueue::buildHeap() {
        for (int i = size / 2 - 1; i >= 0; i--) {
            heapify(i);
        }
    }

}  // namespace Koala
