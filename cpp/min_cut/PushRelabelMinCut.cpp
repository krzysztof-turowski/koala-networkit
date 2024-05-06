#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <algorithm>
#include <tuple>

#include <min_cut/PushRelabelMinCut.hpp> 

namespace Koala {
    
int PushRelabelMinCut::overflowVertex() {
    for (int i = 0; i < numVertices; ++i) {
        if (i != source && i != sink && excess[i] > 0) return i;
    }
    return -1;
}

void PushRelabelMinCut::initializePreflow() {
    height[source] = numVertices;
    for (Edge& e : graph[source]) {
        e.flow = e.capacity;
        graph[e.sink][e.revIndex].flow -= e.capacity;
        excess[e.sink] += e.capacity;
        excess[source] -= e.capacity;
    }
}

void PushRelabelMinCut::push(Edge& e) {
    long long flow = std::min(excess[e.source], e.capacity - e.flow);
    e.flow += flow;
    graph[e.sink][e.revIndex].flow -= flow;
    excess[e.source] -= flow;
    excess[e.sink] += flow;
}

void PushRelabelMinCut::relabel(int v) {
    int minHeight = std::numeric_limits<int>::max();
    for (const Edge& e : graph[v]) {
        if (e.capacity > e.flow) {
            minHeight = std::min(minHeight, height[e.sink]);
        }
    }
    height[v] = minHeight + 1;
}

void PushRelabelMinCut::discharge(int v) {
    while (excess[v] > 0) {
        if (current[v] == graph[v].size()) {
            relabel(v);
            current[v] = 0;
        } else {
            Edge& e = graph[v][current[v]];
            if (e.capacity > e.flow && height[v] > height[e.sink]) {
                push(e);
            }
            current[v]++;
        }
    }
}

void PushRelabelMinCut::solve() {
    initializePreflow();
    int v;
    while ((v = overflowVertex()) != -1) {
        discharge(v);
    }
}

PushRelabelMinCut::PushRelabelMinCut(int vertices, int src, int snk, const std::vector<std::vector<int>>& graphMatrix)
    : numVertices(vertices), source(src), sink(snk), graph(vertices), excess(vertices, 0), height(vertices, 0), active(vertices, false), current(vertices, 0) {
    for (int i = 0; i < vertices; ++i) {
        for (int j = 0; j < vertices; ++j) {
            if (graphMatrix[i][j] > 0) {
                addEdge(i, j, graphMatrix[i][j]);
            }
        }
    }
}

void PushRelabelMinCut::addEdge(int from, int to, long long capacity) {
    Edge forward = {from, to, capacity, 0, static_cast<int>(graph[to].size())};
    Edge reverse = {to, from, 0, 0, static_cast<int>(graph[from].size())};
    graph[from].push_back(forward);
    graph[to].push_back(reverse);
}

long long PushRelabelMinCut::getMaxFlow() const {
    return excess[sink];
}

std::vector<int> PushRelabelMinCut::getMinCut() const {
    std::vector<int> visited(numVertices, 0);
    std::queue<int> queue;

    queue.push(source);
    visited[source] = 1;

    while (!queue.empty()) {
        int v = queue.front();
        queue.pop();

        for (const Edge& e : graph[v]) {
            if (e.flow < e.capacity && !visited[e.sink]) {
                visited[e.sink] = 1;
                queue.push(e.sink);
            }
        }
    }

    return visited;
}

}  // namespace Koala