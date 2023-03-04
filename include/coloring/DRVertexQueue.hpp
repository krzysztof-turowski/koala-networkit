/*
* DRVertexQueue.hpp
*
*  Created on: 28.02.2023
*      Author: Zofia Glapa (zofia.glapa@student.uj.edu.pl)
*/

#pragma once

#include <vector>
#include <utility>
#include <unordered_map>

namespace Koala {

    struct NodeValuePair {
        int node;
        int value;
    };

    bool operator>(const NodeValuePair& a, const NodeValuePair& b);
    bool operator==(const NodeValuePair& a, const NodeValuePair& b);

    class DRVertexQueue {

    public:
        DRVertexQueue();
        DRVertexQueue(const std::vector<std::pair<int, int>> & nodes);
        void updateValue(const int & v, const int& value);
        void insert(const int & v, const int& value);
        NodeValuePair pop();
        bool empty();

    private:
        std::vector<NodeValuePair> queue;
        std::unordered_map<int, int> position;
        int size;

        void swap(const int& i, const int& j);
        void heapify(const int& i);
        void buildHeap();

    };

}